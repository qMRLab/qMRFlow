#!/usr/bin/env nextflow

/*
This workflow contains pre- and post-processing steps to 
calculate longitudinal relaxation time (T1) map using
Magnetization Prepared 2 Rapid Acquisition Gradient Echoes (MP2RAGE)

Dependencies: 
    These dependencies must be installed if Docker is not going
    to be used. 
        - Advanced notmarization tools (ANTs, https://github.com/ANTsX/ANTs)
        - FSL  
        - qMRLab (https://qmrlab.org) 
        - git     

Docker: 
        - https://hub.docker.com/u/qmrlab
        - qmrlab/minimal:v2.3.1
        - qmrlab/antsfsl:latest

Written by: Agah Karakuzu, Juan Jose Velazquez Reyes | 2021
GitHub:     @agahkarakuzu, @jvelazquez-reyes

Users: Please see USAGE for further details
 */

/*Set defaults for parameters determining logic flow to false*/
params.root = false 
params.help = false

/* Call to the mp2rage_wrapper.m will be invoked by params.runcmd.
Depending on the params.platform selection, params.runcmd 
may point to MATLAB or Octave. 
*/
if (params.platform == "octave"){

    if (params.octave_path){
        log.info "Using Octave executable declared in nextflow.config."
        params.octave = params.octave_path + " --no-gui --eval"
    }else{
        log.info "Using Octave in Docker or (if local) from the sys path."
        params.octave = "octave --no-gui --eval"
    }

    params.runcmd = params.octave 
}

if (params.platform == "matlab"){
   
    if (params.matlab_path){
        log.info "Using MATLAB executable declared in nextflow.config."
        params.matlab = params.matlab_path + " -nodisplay -nosplash -nodesktop -r"
    }else{
        log.info "Using MATLAB from the sys path."
        params.matlab = "matlab -nodisplay -nosplash -nodesktop -r"
    }

    params.runcmd = params.matlab
}

params.wrapper_repo = "https://github.com/qMRLab/qMRWrappers.git"
              
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

/*Define bindings for --help*/
if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["use_b1map":"$params.use_b1map",
                "use_bet":"$params.use_bet",
                "bet_recursive":"$params.bet_recursive",
                "bet_threshold":"$params.bet_threshold",
                "platform":"$params.platform",
                "matlab_path":"$params.matlab_path",
                "octave_path":"$params.octave_path",
                "qmrlab_path":"$params.qmrlab_path"
                ]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

/*Scrape file names from a BIDS-compatible dataset
Note:
    BIDS for qMRI is currently under development (BEP001,https://github.com/bids-standard/bep001)
    The current format is valid as of late 2019 and subjected to change.
    For B1plusmaps, (optional) maps are assumed to be located at the fmap
    folder with _TB1map suffix.   
*/
if(params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    
    /* ==== BIDS: MP2RAGE inputs ==== */  
    /* Here, alphabetical indexes matter. */
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_UNIT1.nii.gz", maxDepth: 2, size: 1, flat: true)
    (unit) = in_data
        .map{sid, unit1  -> [tuple(sid, unit1)]}                                   
        .separate(1)

    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_UNIT1.json", maxDepth: 2, size: 1, flat: true)
    (unitj) = in_data
        .map{sid, unit1  -> [tuple(sid, unit1)]}                                   
        .separate(1)

    /* ==== BIDS: B1 map ==== */             
    /* Look for B1map in fmap folder */
    b1_data = Channel
           .fromFilePairs("$root/**/fmap/sub-*_TB1map.nii.gz", maxDepth:2, size:1, flat:true)   
    (b1raw) = b1_data       
           .map{sid, B1plusmap -> [tuple(sid, B1plusmap)]}     
           .separate(1)  
}   
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 

log.info "qMRflow: MP2RAGE pipeline"
log.info "======================="
log.info ""
log.info "##     ## #######   ######  #######     ###      #####  ########"
log.info "###   ### ##    ## ###  ### ##    ##   ## ##    ##   ## ##"
log.info "#### #### ##    ## ##  ###  ##    ##  ##   ##  ##    ## ##"
log.info "## ### ## #######     ###   #######  ##     ## ##       #####"
log.info "##     ## ##         ###    ##  ##   ######### ##  #### ##"
log.info "##     ## ##        ###     ##   ##  ##     ##  ##   ## ##"
log.info "##     ## ##       #######  ##    ## ##     ##   #####  ########"
log.info ""
log.info "Start time: $workflow.start"
log.info ""
log.info ""
log.info "DATA"
log.info "===="
log.info ""
log.info "BIDS option has been enabled."
log.warn "qMRI protocols will be read from sidecar .json files."
log.info ""
log.info "OPTIONS"
log.info "======="
log.info ""
log.info "[GLOBAL]"
log.info "---------------"
log.info "Selected platform: $params.platform"
log.info "BET enabled: $params.use_bet"
log.info "B1+ correction enabled: $params.use_b1map"
log.info ""
log.info "[FSL BET]"
log.info "---------------"
log.info "Enabled: $params.use_bet"
log.info "Fractional intensity threshold: $params.bet_threshold"
log.info "Robust brain center estimation: $params.bet_recursive"
log.info ""
log.info "[qMRLab mp2rage]"
log.info "---------------"
log.warn "Acquisition protocols will be read from  sidecar .json files (BIDS)."
if (params.use_b1map){
log.info "B1+ correction has been ENABLED."  
log.warn "Process will be skipped for participants missing a B1map file."}
if (!params.use_b1map){
log.info "B1+ correction has been DISABLED."
log.warn "Process will NOT take any (possibly) existing B1maps into account."
}
log.info ""
log.info "======================="

/* Split input data into two to deal with BET and B1map cases */
unit.into{unit_ch1;unit_ch2;unit_ch3}

unitj.into{unitj_ch1;unitj_ch2}

process Extract_Brain{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(unit) from unit_ch3

    output:
        tuple val(sid), "${sid}_UNIT1_mask.nii.gz" optional true into mask_from_bet
        file "${sid}_UNIT1_mask.nii.gz"

    script:
         if (params.bet_recursive){
        """    
        bet $unit ${sid}_UNIT1.nii.gz -m -R -n -f $params.bet_threshold
        """}
        else{
        """    
        bet $unit ${sid}_UNIT1.nii.gz -m -n -f $params.bet_threshold
        """
        }

}

/* There is no input optional true concept in nextflow
The process consuming the individual input channels will 
only execute if the channel is populated.
*/

/* We need to conditionally create channels with
input data or as empty channels.
*/

if (!params.use_bet){
    Channel
        .empty()
        .set{mask_from_bet}
}

/* Split mask_from_bet into two to deal with B1map cases later. */
mask_from_bet.into{mask_from_bet_ch1;mask_from_bet_ch2}

unit_ch1
    .join(unitj_ch1)
    .set{mp2rage_without_b1}

/* Split into two to deal with mask cases. */
mp2rage_without_b1.into{mp2rage_wo_b1_bet;mp2rage_wo_b1}

unit_ch2
    .join(unitj_ch2)
    .join(b1raw)
    .set{mp2rage_with_b1}

/* Split into two to deal with mask cases. */
mp2rage_with_b1.into{mp2rage_w_b1_bet;mp2rage_w_b1}

/* Deal with mask cases. */
mp2rage_wo_b1_bet
    .join(mask_from_bet_ch1)
    .set{mp2rage_wo_b1_bet_merged}

mp2rage_w_b1_bet
    .join(mask_from_bet_ch2)
    .set{mp2rage_w_b1_bet_merged}

/* Depending on the nextflow.config 
settings for use_b1map and use_bet, one of the
following 4 processes will be executed. 
*/

process Fit_MP2RAGE_With_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_b1map == true && params.use_bet == true

    input:
        tuple val(sid), file(uni), file(unij), file(b1map), file(mask) from mp2rage_w_b1_bet_merged
        
    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_T1map.json"
	file "${sid}_UNIT1.nii.gz"
	file "${sid}_UNIT1.json"
        file "${sid}_mp2rage.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            git clone -b mp2rage_UNIT1 $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mp2rage_UNIT1_wrapper('$uni','$unij','mask','$mask','b1map','$b1map','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

        """
}

process Fit_MP2RAGE_With_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_b1map == true && params.use_bet == false

    input:
        tuple val(sid), file(uni), file(unij), file(b1map) from mp2rage_w_b1
        
    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_T1map.json"
	file "${sid}_UNIT1.nii.gz"
	file "${sid}_UNIT1.json"
        file "${sid}_mp2rage.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            git clone -b mp2rage_UNIT1 $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mp2rage_UNIT1_wrapper('$uni','$unij','b1map','$b1map','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

        """
}

process Fit_MP2RAGE_Without_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_b1map == false && params.use_bet == true

    input:
        tuple val(sid), file(uni), file(unij), file(mask) from mp2rage_wo_b1_bet_merged
        
    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_T1map.json"
	file "${sid}_UNIT1.nii.gz"
	file "${sid}_UNIT1.json"
        file "${sid}_mp2rage.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            git clone -b mp2rage_UNIT1 $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mp2rage_UNIT1_wrapper('$uni','$unij','mask','$mask','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

        """
}

process Fit_MP2RAGE_Without_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'
    
    when:
        params.use_b1map == false && params.use_bet == false

    input:
        tuple val(sid), file(uni), file(unij) from mp2rage_wo_b1

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_UNIT1.nii.gz" 
        file "${sid}_T1map.json" 
        file "${sid}_UNIT1.json" 
        file "${sid}_mp2rage.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            git clone -b mp2rage_UNIT1 $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mp2rage_UNIT1_wrapper('$uni','$unij','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

        """
}

