#!/usr/bin/env nextflow

/*
This workflow contains pre- and post-processing steps to 
calculate Actual Flip-Angle Imaging (AFI) b1+ map.

Dependencies: 
    These dependencies must be installed if Docker is not going
    to be used. 
        - Advanced notmarization tools (ANTs, https://github.com/ANTsX/ANTs)
        - FSL  
        - qMRLab (https://qmrlab.org) > 2.4.1
        - git     

Docker: 
        - AFI model is not defined yet in version 2.4.1

Author:
    Agah Karakuzu and Juan José Velázquez Reyes 2021
    agahkarakuzu@gmail.com 

Users: Please see USAGE for further details
 */

/*Set defaults for parameters determining logic flow to false*/
params.root = false 
params.help = false

/* Call to the b1_afi_wrapper.m will be invoked by params.runcmd.
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
    bindings = ["ants_dim":"$params.ants_dim",
                "ants_metric":"$params.ants_metric",
                "ants_metric_weight":"$params.ants_metric_weight",
                "ants_metric_bins":"$params.ants_metric_bins",
                "ants_metric_sampling":"$params.ants_metric_sampling",
                "ants_metric_samplingprct":"$params.ants_metric_samplingprct",
                "ants_transform":"$params.ants_transform",
                "ants_convergence":"$params.ants_convergence",
                "ants_shrink":"$params.ants_shrink",
                "ants_smoothing":"$params.ants_smoothing",
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
    For B1plusmaps, there is not a specification yet. To circumvent this 
    issue, these (optional) maps are assumed to be located at the fmap
    folder with _TB1map suffix.   
*/
if(params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    
    /* ==== BIDS: b1afi inputs ==== */  
    /* Here, alphabetical indexes matter. Therefore, tr1 -> tr2 */
    in_data = Channel
        .fromFilePairs("$root/**/fmap/sub-*_acq-tr{1,2}_TB1AFI.nii.gz", maxDepth: 3, size: 2, flat: true)
    (afiData1, afiData2) = in_data
        .map{sid, tr1, tr2  -> [    tuple(sid, tr1),
                                    tuple(sid, tr2)]}                                   
        .separate(2)

    in_data = Channel
        .fromFilePairs("$root/**/fmap/sub-*_acq-tr{1,2}_TB1AFI.json", maxDepth: 3, size: 2, flat: true)
    (afiData1j, afiData2j) = in_data
        .map{sid, tr1, tr2  -> [    tuple(sid, tr1),
                                    tuple(sid, tr2)]}                                   
        .separate(2)    
}   
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 

/*Split T1w into three channels
    afiData1_pre_ch2 --> afiData1_for_bet
    afiData1_pre_ch3 --> afiData1_post
*/
afiData1.into{afiData1_for_bet; afiData1_post}

log.info "qMRflow: b1afi pipeline"
log.info "======================="
log.info ""
log.info "#######    ##        ###    ######## ##"
log.info "##    ##  ###       ## ##   ##       ##"
log.info "##    ## ####      ##   ##  ##       ##"
log.info "#######    ## ### ##     ## #####    ##"
log.info "##    ##   ##     ######### ##       ##"
log.info "##    ##   ##     ##     ## ##       ##"
log.info "#######    ##     ##     ## ##       ##"
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
log.info ""
log.info "[ANTs Registration]"
log.info "-------------------"
log.info "Dimensionality: $params.ants_dim"
log.info "Metric: $params.ants_metric"
log.info "Weight: $params.ants_metric_weight"
log.info "Number of bins: $params.ants_metric_bins"
log.info "Sampling type: $params.ants_metric_sampling"
log.info "Sampling percentage: $params.ants_metric_samplingprct"
log.info "Transform: $params.ants_transform"
log.info "Convergence: $params.ants_convergence"
log.info "Shrink factors: $params.ants_shrink"
log.info "Smoothing sigmas: $params.ants_smoothing"
log.info ""
log.info "[FSL BET]"
log.info "---------------"
log.info "Enabled: $params.use_bet"
log.info "Fractional intensity threshold: $params.bet_threshold"
log.info "Robust brain center estimation: $params.bet_recursive"
log.info ""
log.info "[qMRLab b1_afi]"
log.info "---------------"
log.warn "Acquisition protocols will be read from  sidecar .json files (BIDS)."
log.info ""
log.info "======================="

process Extract_Brain{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(afiData1) from afiData1_for_bet

    output:
        tuple val(sid), "${sid}_acq-afiData1_mask.nii.gz" optional true into mask_from_bet
        file "${sid}_acq-afiData1_mask.nii.gz"

    script:
         if (params.bet_recursive){
        """    
        bet $afiData1 ${sid}_acq-afiData1.nii.gz -m -R -n -f $params.bet_threshold
        """}
        else{
        """    
        bet $afiData1 ${sid}_acq-afiData1.nii.gz -m -n -f $params.bet_threshold
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

/* Split mask_from_bet into two for two Fiting processes. */
mask_from_bet.into{mask_from_bet_ch1;mask_from_bet_ch2}

/*Merge afiData1_post with afiData2.*/
/*Without mask */
afiData1_post
    .join(afiData2)
    .join(afiData1j)
    .join(afiData2j)
    .set{b1afi_for_fitting}

b1afi_for_fitting.into{b1afi_without_mask;b1afi_for_fitting_ch2}

/* With mask */
b1afi_for_fitting_ch2
    .join(mask_from_bet_ch1)
    .set{b1afi_with_mask}


/* Depending on the nextflow.config 
settings for use_bet, one of the
following 2 processes will be executed. 
*/

process Fit_b1afi_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*TB1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(afiData1), file(afiData2), file(afiData1j), file(afiData2j), file(mask) from b1afi_with_mask
        
    output:
	tuple val(sid), "${sid}_TB1map.nii.gz"
        file "${sid}_TB1map.nii.gz" 
        file "${sid}_TB1map.json" 
        file "${sid}_b1_afi.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/b1_afi/b1_afi_wrapper.m b1_afi_wrapper.m

            $params.runcmd "b1_afi_wrapper('$afiData1','$afiData2','$afiData1j','$afiData2j','mask','$mask','b1filter','$params.b1_filter','type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path','sid','${sid}','containerType','$workflow.containerEngine', 'containerTag','$params.containerTag','description','$params.description','datasetDOI','$params.datasetDOI','datasetURL','$params.datasetURL','datasetVersion','$params.datasetVersion'); exit();"

        """
}

process Fit_b1afi_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.nii.gz'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*TB1map.json'
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy', pattern: '*.mat'
    publishDir "$root/derivatives/qMRLab", mode: 'copy', pattern: 'dataset_description.json'

    when:
        params.use_bet == false

    input:
        tuple val(sid), file(afiData1), file(afiData2), file(afiData1j), file(afiData2j) from b1afi_without_mask

    output:
	tuple val(sid), "${sid}_TB1map.nii.gz"
        file "${sid}_TB1map.nii.gz"  
        file "${sid}_TB1map.json" 
        file "${sid}_b1_afi.qmrlab.mat"
	file "dataset_description.json"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/b1_afi/b1_afi_wrapper.m b1_afi_wrapper.m

            $params.runcmd "b1_afi_wrapper('$afiData1','$afiData2','$afiData1j','$afiData2j','b1filter','$params.b1_filter','type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size, 'qmrlab_path','$params.qmrlab_path','sid','${sid}','containerType','$workflow.containerEngine', 'containerTag','$params.containerTag','description','$params.description','datasetDOI','$params.datasetDOI','datasetURL','$params.datasetURL','datasetVersion','$params.datasetVersion'); exit();"

        """
               
}

