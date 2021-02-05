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
    (uni) = in_data
        .map{sid, unit1  -> [tuple(sid, unit1)]}                                   
        .separate(1)

    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_UNIT1.json", maxDepth: 2, size: 1, flat: true)
    (unij) = in_data
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

/* Split input data into two to deal with B1map cases */
uni.into{uni_ch1;uni_ch2}

unij.into{unij_ch1;unij_ch2}

uni_ch1
    .join(unij_ch1)
    .set{mp2rage_without_b1}

uni_ch2
    .join(unij_ch2)
    .join(b1raw)
    .set{mp2rage_with_b1}

/* Depending on the nextflow.config 
settings for use_b1map, one of the
following 2 processes will be executed. 
*/

process Fit_MP2RAGE_With_B1map{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_b1map == true

    input:
        tuple val(sid), file(uni), file(unij), file(b1map) from mp2rage_with_b1
        
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

process Fit_MP2RAGE_Without_B1map{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*T1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*UNIT1.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'
    
    when:
        params.use_b1map == false

    input:
        tuple val(sid), file(uni), file(unij) from mp2rage_without_b1

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

