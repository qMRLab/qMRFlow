#!/usr/bin/env nextflow

/*
This workflow contains pre- and post-processing steps to 
calculate Magnetization Transfer Saturation Index (MTsat) map along
with a longitudinal relaxation time (T1) map.

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

Author:
    Agah Karakuzu 2019
    agahkarakuzu@gmail.com 

Users: Please see USAGE for further details
 */

/*Set defaults for parameters determining logic flow to false*/
params.root = false 
params.help = false

/* Call to the mt_sat_wrapper.m will be invoked by params.runcmd.
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
    folder with _B1plusmap suffix.   
*/
if(params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    
    /* ==== BIDS: MTRinputs ==== */  
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_acq-{MToff,MTon}_MTR.nii*", maxDepth: 2, size: 2, flat: true)
    (mtoff, mton) = in_data
        .map{sid, MToff, MTon-> [    tuple(sid, MToff),
                                            tuple(sid, MTon)]}                                   
        .separate(2)

}   
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

mton.into{mton_for_bet; mton_ch1;mton_ch2; mton_ch3}

/* Pair for alignment */ 
mton_ch1
    .join(mtoff)
    .set{mtr_for_alignment}

log.info "qMRflow: MTR pipeline"
log.info "======================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""
log.info ""
log.info "DATA"
log.info "===="
log.info ""
log.info "BIDS option has been enabled."
log.warn "No protocol will be read for this operaiton."
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
log.info "[qMRLab mt_ratio]"
log.info "---------------"
log.info "No options or protocol are available for this model."
log.info "https://qmrlab.readthedocs.io/en/master/mt_ratio_batch.html"
log.info ""
log.info "======================="

process Align_Input_Volumes {
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    input:
        tuple val(sid), file(mton), file(mtoff) from mtr_for_alignment

    output:
        tuple val(sid), "${sid}_acq-MToff_MTR_aligned.nii.gz" into mtoff_aligned
        file "${sid}_acq-MToff_MTR_aligned.nii.gz"
        file "${sid}_mtoff_to_mton_displacement.*.mat"

    script:
        """
        antsRegistration -d $params.ants_dim \
                            --float 0 \
                            -o [${sid}_mtoff_to_mton_displacement.mat,${sid}_acq-MToff_MTR_aligned.nii.gz] \
                            --transform $params.ants_transform \
                            --metric $params.ants_metric[$mton,$mtoff,$params.ants_metric_weight, $params.ants_metric_bins,$params.ants_metric_sampling,$params.ants_metric_samplingprct] \
                            --convergence $params.ants_convergence \
                            --shrink-factors $params.ants_shrink \
                            --smoothing-sigmas $params.ants_smoothing
        """
}


process Extract_Brain{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(mton) from mton_for_bet

    output:
        tuple val(sid), "${sid}_acq-MTon_mask.nii.gz" optional true into mask_from_bet
        file "${sid}_acq-MTon_mask.nii.gz"

    script:
         if (params.bet_recursive){
        """    
        bet $mton ${sid}_acq-MTon.nii.gz -m -R -n -f $params.bet_threshold
        """}
        else{
        """    
        bet $mton ${sid}_acq-MTon.nii.gz -m -n -f $params.bet_threshold
        """
        }

}

if (!params.use_bet){
    Channel
        .empty()
        .set{mask_from_bet}
}

mtoff_aligned.into{mtoff_aligned_ch1;mtoff_aligned_ch2}

mton_ch2
    .join(mtoff_aligned_ch1)
    .join(mask_from_bet)
    .set{mtr_with_bet}

process Fit_MTR_With_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(mton), file(mtoff), file(mask) from mtr_with_bet
        
    output:
        file "${sid}_MTRmap.nii.gz" 
        file "${sid}_MTRmap.json" 
        file "${sid}_mt_ratio.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_ratio_wrapper('$mton','$mtoff','mask','$mask','qmrlab_path','$params.qmrlab_path', 'sid','${sid}'); exit();"
        """
}

mton_ch3
    .join(mtoff_aligned_ch2)
    .set{mtr_without_bet}

process Fit_MTR_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == false

    input:
        tuple val(sid), file(mton), file(mtoff) from mtr_without_bet
        
    output:
        file "${sid}_MTRmap.nii.gz" 
        file "${sid}_MTRmap.json" 
        file "${sid}_mt_ratio.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_ratio_wrapper('$mton','$mtoff','qmrlab_path','$params.qmrlab_path', 'sid','${sid}'); exit();"
        """
}