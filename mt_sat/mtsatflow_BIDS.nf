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
                "use_b1cor":"$params.use_b1cor",
                "b1cor_factor":"$params.b1cor_factor",
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
    
    /* ==== BIDS: MTSat inputs ==== */  
    /* Here, alphabetical indexes matter. Therefore, MToff -> MTon -> T1w */
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_flip-{01,02}_mt-{off,on}_MTS.nii.gz", maxDepth: 3, size: 3, flat: true)
    (pdw, mtw, t1w) = in_data
        .map{sid, MToff, MTon, T1w  -> [    tuple(sid, MToff),
                                            tuple(sid, MTon),
                                            tuple(sid, T1w)]}                                   
        .separate(3)

    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_flip-{01,02}_mt-{off,on}_MTS.json", maxDepth: 3, size: 3, flat: true)
    (pdwj, mtwj, t1wj) = in_data
        .map{sid, MToff, MTon, T1w  -> [    tuple(sid, MToff),
                                            tuple(sid, MTon),
                                            tuple(sid, T1w)]}                                   
        .separate(3)    

    /* ==== BIDS: B1 map ==== */             
    /* Look for B1map in fmap folder */
    b1_data = Channel
           .fromFilePairs("$root/**/fmap/sub-*_TB1map.nii.gz", maxDepth:2, size:1, flat:true)   
    (b1raw) = b1_data       
           .map{sid, TB1map -> [tuple(sid, TB1map)]}     
           .separate(1)  
}   
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

/*Each data type is defined as a channel. To pass all the channels 
  to the same process accurately, these channels must be joined. 
*/ 

/*Split T1w into three channels
    t1w_pre_ch1 --> mtsat_for_alignment
    t1w_pre_ch2 --> t1w_for_bet
    t1w_pre_ch3 --> t1w_post
*/
t1w.into{t1w_pre_ch1; t1w_for_bet; t1w_post}

/* Merge PDw, MTw and T1w for alignment*/
pdw 
    .join(mtw)
    .join(t1w_pre_ch1)
    .set{mtsat_for_alignment}

log.info "qMRflow: MTsat pipeline"
log.info "======================="
log.info ""
log.info "##     ## ########  ######     ###    ########"
log.info "###   ###    ##    ##    ##   ## ##      ##"
log.info "#### ####    ##    ##        ##   ##     ##"
log.info "## ### ##    ##     ######  ##     ##    ##"
log.info "##     ##    ##          ## #########    ##"
log.info "##     ##    ##    ##    ## ##     ##    ##"
log.info "##     ##    ##     ######  ##     ##    ##"
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
log.info "B1+ correction enabled: $params.use_b1cor"
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
log.info "[qMRLab mt_sat]"
log.info "---------------"
log.warn "Acquisition protocols will be read from  sidecar .json files (BIDS)."
if (params.use_b1cor){
log.info "B1+ correction has been ENABLED."  
log.warn "Process will be skipped for participants missing a B1map file."   
log.info "B1 correction factor: $params.b1cor_factor"}
if (!params.use_b1cor){
log.info "B1+ correction has been DISABLED."
log.warn "Process will NOT take any (possibly) existing B1maps into account."
}
log.info ""
log.info "======================="

/*Perform rigid registration to correct for head movement across scans:
    - MTw (moving) --> T1w (fixed)
    - PDw (moving) --> T1w (fixed)
*/     

process Align_Input_Volumes {
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    input:
        tuple val(sid), file(pdw), file(mtw), file(t1w) from mtsat_for_alignment

    output:
        tuple val(sid), "${sid}_acq-MTon_MTS_aligned.nii.gz", "${sid}_acq-MToff_MTS_aligned.nii.gz"\
        into mtsat_from_alignment
        file "${sid}_acq-MTon_MTS_aligned.nii.gz"
        file "${sid}_acq-MToff_MTS_aligned.nii.gz"
        file "${sid}_mtw_to_t1w_displacement.*.mat"
        file "${sid}_pdw_to_t1w_displacement.*.mat"

    script:
        """
        antsRegistration -d $params.ants_dim \
                            --float 0 \
                            -o [${sid}_mtw_to_t1w_displacement.mat,${sid}_acq-MTon_MTS_aligned.nii.gz] \
                            --transform $params.ants_transform \
                            --metric $params.ants_metric[$t1w,$mtw,$params.ants_metric_weight, $params.ants_metric_bins,$params.ants_metric_sampling,$params.ants_metric_samplingprct] \
                            --convergence $params.ants_convergence \
                            --shrink-factors $params.ants_shrink \
                            --smoothing-sigmas $params.ants_smoothing

        antsRegistration -d $params.ants_dim \
                            --float 0 \
                            -o [${sid}_pdw_to_t1w_displacement.mat,${sid}_acq-MToff_MTS_aligned.nii.gz] \
                            --transform $params.ants_transform \
                            --metric $params.ants_metric[$t1w,$pdw,$params.ants_metric_weight, $params.ants_metric_bins,$params.ants_metric_sampling,$params.ants_metric_samplingprct] \
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
        tuple val(sid), file(t1w) from t1w_for_bet

    output:
        tuple val(sid), "${sid}_acq-T1w_mask.nii.gz" optional true into mask_from_bet
        file "${sid}_acq-T1w_mask.nii.gz"

    script:
         if (params.bet_recursive){
        """    
        bet $t1w ${sid}_acq-T1w.nii.gz -m -R -n -f $params.bet_threshold
        """}
        else{
        """    
        bet $t1w ${sid}_acq-T1w.nii.gz -m -n -f $params.bet_threshold
        """
        }

}

/* Split t1w_post into two to deal with B1map cases */
t1w_post.into{t1w_post_ch1;t1w_post_ch2; t1w_post_ch3}

/* Split mtsat_from_alignment into two to deal with B1map cases */
mtsat_from_alignment.into{mfa_ch1;mfa_ch2}

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
mask_from_bet.into{mask_from_bet_ch1;mask_from_bet_ch2;mask_from_bet_ch3}

t1wj.into{t1wj_ch1;t1wj_ch2}
pdwj.into{pdwj_ch1;pdwj_ch2}
mtwj.into{mtwj_ch1;mtwj_ch2}

t1w_post_ch3
    .join(b1raw)
    .set{b1_for_alignment}

/* For now, just use identity transform (i.e. upsample w/o additional transformation).*/

process B1_Align{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true

    input:
        tuple val(sid), file(t1w), file(b1raw) from b1_for_alignment
        

    output:
        tuple val(sid), "${sid}_B1plusmap_aligned.nii.gz" optional true into b1_aligned
        file "${sid}_B1plusmap_aligned.nii.gz"

    script:
        """
        antsApplyTransforms -d 3 -e 0 -i $b1raw \
                            -r $t1w \
                            -o ${sid}_B1plusmap_aligned.nii.gz \
                            -t identity
        """

}

if (!params.use_b1cor){
    Channel
        .empty()
        .set{b1_aligned}
}

b1_aligned.into{b1_aligned_ch1;b1_for_smooth_without_mask}

b1_aligned_ch1
   .join(mask_from_bet_ch3)
   .set{b1_for_smooth_with_mask}
            
process B1_Smooth_With_Mask{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    if (!params.matlab_path_exception){
    container 'qmrlab/minimal:v2.3.1'
    }
    
    when:
        params.use_b1cor == true && params.use_bet == true

    input:
        tuple val(sid), file(b1aligned), file(mask) from b1_for_smooth_with_mask

    output:
        tuple val(sid), "${sid}_B1plusmap_filtered.nii.gz" optional true into b1_filtered_w_mask 
        file "${sid}_B1plusmap_filtered.nii.gz"
        file "${sid}_B1plusmap_filtered.json"

    script: 
        if (params.matlab_path_exception){
        """
            git clone -b mt_sat-argparser $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.matlab_path_exception -nodesktop -nosplash -r "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned', 'mask', '$mask', 'type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path_exception','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """
        }else{
        """
            git clone -b mt_sat-argparser $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned', 'mask', '$mask', 'type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """

        }

}

process B1_Smooth_Without_Mask{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    container 'qmrlab/minimal:v2.3.1'

    when:
        params.use_b1cor == true && params.use_bet == false

    input:
        tuple val(sid), file(b1aligned) from b1_for_smooth_without_mask
    
    output:
        tuple val(sid), "${sid}_B1plusmap_filtered.nii.gz" optional true into b1_filtered_wo_mask 
        file "${sid}_B1plusmap_filtered.nii.gz"
        file "${sid}_B1plusmap_filtered.json"
        
    script:
    if (params.matlab_path_exception){
        """
            git clone -b mt_sat-argparser $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.matlab_path_exception -nodesktop -nosplash -r "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned','type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path_exception','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """
        }else{
        """
            git clone -b mt_sat-argparser $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..
            
            $params.runcmd "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned', 'type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """

    }

}

if (params.use_bet){
/*Merge tw1_post with mtsat_from_alignment and b1plus.*/
t1w_post_ch1
    .join(mfa_ch1)
    .join(b1_filtered_w_mask)
    .join(t1wj_ch1)
    .join(mtwj_ch1)
    .join(pdwj_ch1)
    .set{mtsat_for_fitting_with_b1}
}else{
t1w_post_ch1
    .join(mfa_ch1)
    .join(b1_filtered_wo_mask)
    .join(t1wj_ch1)
    .join(mtwj_ch1)
    .join(pdwj_ch1)
    .set{mtsat_for_fitting_with_b1}

}



mtsat_for_fitting_with_b1.into{mtsat_with_b1_bet;mtsat_with_b1}

/*Merge tw1_post with mtsat_from_alignment only.
WITHOUT B1 MAP
*/
t1w_post_ch2
    .join(mfa_ch2)
    .join(t1wj_ch2)
    .join(mtwj_ch2)
    .join(pdwj_ch2)
    .set{mtsat_for_fitting_without_b1}

mtsat_for_fitting_without_b1.into{mtsat_without_b1_bet;mtsat_without_b1}

/* We need to join these channels to avoid problems.
WITH B1 MAP
*/
mtsat_with_b1_bet
    .join(mask_from_bet_ch1)
    .set{mtsat_with_b1_bet_merged}

/* Depeding on the nextflow.config 
settings for use_b1cor and use_bet, one of th
following 4 processes will be executed. 
*/

process Fit_MTsat_With_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true && params.use_bet == true

    input:
        tuple val(sid), file(t1w), file(mtw_reg), file(pdw_reg),\
        file(b1map), file(t1wj), file(mtwj), file(pdwj), file(mask) from mtsat_with_b1_bet_merged
        
    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_MTsat.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_MTsat.json"  
        file "${sid}_mt_sat.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_sat_wrapper('$mtw_reg','$pdw_reg','$t1w','$mtwj','$pdwj','$t1wj','mask','$mask','b1map','$b1map','b1factor',$params.b1cor_factor,'qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}

process Fit_MTsat_With_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true && params.use_bet == false

    input:
        tuple val(sid), file(t1w), file(mtw_reg), file(pdw_reg),\
        file(b1map), file(t1wj), file(mtwj), file(pdwj) from mtsat_with_b1

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_MTsat.nii.gz" 
        file "${sid}_T1map.json" 
        file "${sid}_MTsat.json" 
        file "${sid}_mt_sat.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_sat_wrapper('$mtw_reg','$pdw_reg','$t1w','$mtwj','$pdwj','$t1wj','b1map','$b1map','b1factor',$params.b1cor_factor,'qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
               
}


/* We need to join these channels to avoid problems.*/
mtsat_without_b1_bet
    .join(mask_from_bet_ch2)
    .set{mtsat_without_b1_bet_merged}

process Fit_MTsat_Without_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == false && params.use_bet==true

    input:
        tuple val(sid), file(t1w), file(mtw_reg), file(pdw_reg),\
        file(t1wj), file(mtwj), file(pdwj), file(mask) from mtsat_without_b1_bet_merged

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_MTsat.nii.gz" 
        file "${sid}_T1map.json" 
        file "${sid}_MTsat.json" 
        file "${sid}_mt_sat.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_sat_wrapper('$mtw_reg','$pdw_reg','$t1w','$mtwj','$pdwj','$t1wj','mask','$mask','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}

process Fit_MTsat_Without_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == false && params.use_bet==false

    input:
        tuple val(sid), file(t1w), file(mtw_reg), file(pdw_reg),\
        file(t1wj), file(mtwj), file(pdwj) from mtsat_without_b1

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_MTsat.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_MTsat.json"  
        file "${sid}_mt_sat.qmrlab.mat"

    script: 
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.runcmd "addpath(genpath('qMRWrappers')); mt_sat_wrapper('$mtw_reg','$pdw_reg','$t1w','$mtwj','$pdwj','$t1wj','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}
