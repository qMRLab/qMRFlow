#!/usr/bin/env nextflow

/*
This workflow contains pre- and post-processing steps to 
calculate Variable Flip Angle T1 (VFAT1) map

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

/* Call to the vfa_t1_wrapper.m will be invoked by params.runcmd.
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
    
    /* ==== BIDS: VFAT1 inputs ==== */  
    /* Here, alphabetical indexes matter. Therefore, MToff -> MTon -> T1w */
    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_fa-{1,2}_VFA.nii.gz", maxDepth: 4, size: 2, flat: true)
    (vfa1, vfa2) = in_data
        .map{sid, vfa1, vfa2  -> [    tuple(sid, vfa1),
                                      tuple(sid, vfa2)]}                                   
        .separate(2)

    in_data = Channel
        .fromFilePairs("$root/**/anat/sub-*_fa-{1,2}_VFA.json", maxDepth: 4, size: 2, flat: true)
    (vfa1j, vfa2j) = in_data
        .map{sid, vfa1, vfa2  -> [    tuple(sid, vfa1),
                                      tuple(sid, vfa2)]}                                   
        .separate(2)    

    /* ==== BIDS: B1 map ==== */             
    /* Look for B1map in fmap folder */
    b1_data = Channel
           .fromFilePairs("$root/**/fmap/sub-*_TB1map.nii.gz", maxDepth:4, size:1, flat:true)   
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

/*Split vfa1 into three channels
    vfa1_pre_ch1 --> vfat1_for_alignment
    vfa1_pre_ch2 --> vfa1_for_bet
    vfa1_pre_ch3 --> vfa1_post
*/
vfa1.into{vfa1_pre_ch1; vfa1_for_bet; vfa1_post}

/* Merge PDw, MTw and T1w for alignment*/
vfa2 
    .join(vfa1_pre_ch1)
    .set{vfat1_for_alignment}

log.info "qMRflow: VFAT1 pipeline"
log.info "======================="
log.info ""
log.info "##     ## ########   ###   ########  ##"
log.info "##     ## ##        ## ##     ##    ###"
log.info "###   ### ##       ##   ##    ##     ##"
log.info " ##   ##  ####    ##     ##   ##     ##"
log.info " ### ###  ##      #########   ##     ##"
log.info "  #####   ##      ##     ##   ##     ##"
log.info "   ###    ##      ##     ##   ##     ##"
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
log.info "[qMRLab vfa_t1]"
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
    - vfa2 (moving) --> vfa1 (fixed)
*/     

process Align_Input_Volumes {
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    input:
        tuple val(sid), file(vfa2), file(vfa1) from vfat1_for_alignment

    output:
        tuple val(sid), "${sid}_fa-2_VFA_aligned.nii.gz"\
        into vfat1_from_alignment
        file "${sid}_fa-2_VFA_aligned.nii.gz"
        file "${sid}_vfa2_to_vfa1_displacement.*.mat"

    script:
        """
        antsRegistration -d $params.ants_dim \
                            --float 0 \
                            -o [${sid}_vfa2_to_vfa1_displacement.mat,${sid}_fa-2_VFA_aligned.nii.gz] \
                            --transform $params.ants_transform \
                            --metric $params.ants_metric[$vfa1,$vfa2,$params.ants_metric_weight, $params.ants_metric_bins,$params.ants_metric_sampling,$params.ants_metric_samplingprct] \
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
        tuple val(sid), file(vfa1) from vfa1_for_bet

    output:
        tuple val(sid), "${sid}_fa-1_mask.nii.gz" optional true into mask_from_bet
        file "${sid}_fa-1_mask.nii.gz"

    script:
         if (params.bet_recursive){
        """    
        bet $vfa1 ${sid}_fa-1.nii.gz -m -R -n -f $params.bet_threshold
        """}
        else{
        """    
        bet $vfa1 ${sid}_fa-1.nii.gz -m -n -f $params.bet_threshold
        """
        }

}

/* Split vfa1_post into two to deal with B1map cases */
vfa1_post.into{vfa1_post_ch1;vfa1_post_ch2; vfa1_post_ch3}

/* Split vfat1_from_alignment into two to deal with B1map cases */
vfat1_from_alignment.into{vfa_ch1;vfa_ch2}

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

vfa1j.into{vfa1j_ch1;vfa1j_ch2}
vfa2j.into{vfa2j_ch1;vfa2j_ch2}

vfa1_post_ch3
    .join(b1raw)
    .set{b1_for_alignment}

/* For now, just use identity transform (i.e. upsample w/o additional transformation).*/

process B1_Align{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true

    input:
        tuple val(sid), file(vfa1), file(b1raw) from b1_for_alignment
        

    output:
        tuple val(sid), "${sid}_TB1map_aligned.nii.gz" optional true into b1_aligned
        file "${sid}_TB1map_aligned.nii.gz"

    script:
        """
        antsApplyTransforms -d 3 -e 0 -i $b1raw \
                            -r $vfa1 \
                            -o ${sid}_TB1map_aligned.nii.gz \
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
        tuple val(sid), "${sid}_TB1map_filtered.nii.gz" optional true into b1_filtered_w_mask 
        file "${sid}_TB1map_filtered.nii.gz"
        file "${sid}_TB1map_filtered.json"

    script: 
        if (params.matlab_path_exception){
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.matlab_path_exception -nodesktop -nosplash -r "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned', 'mask', '$mask', 'type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path_exception','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """
        }else{
        """
            git clone $params.wrapper_repo 
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
        tuple val(sid), "${sid}_TB1map_filtered.nii.gz" optional true into b1_filtered_wo_mask 
        file "${sid}_TB1map_filtered.nii.gz"
        file "${sid}_TB1map_filtered.json"
        
    script:
    if (params.matlab_path_exception){
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..

            $params.matlab_path_exception -nodesktop -nosplash -r "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned','type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path_exception','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """
        }else{
        """
            git clone $params.wrapper_repo 
            cd qMRWrappers
            sh init_qmrlab_wrapper.sh $params.wrapper_version 
            cd ..
            
            $params.runcmd "addpath(genpath('qMRWrappers')); filter_map_wrapper('$b1aligned', 'type','$params.b1_filter_type','order',$params.b1_filter_order,'dimension','$params.b1_filter_dimension','size',$params.b1_filter_size,'qmrlab_path','$params.qmrlab_path','siemens','$params.b1_filter_siemens', 'sid','${sid}'); exit();" 
        """

    }

}

if (params.use_bet){
/*Merge tw1_post with vfat1_from_alignment and b1plus.*/
vfa1_post_ch1
    .join(vfa_ch1)
    .join(b1_filtered_w_mask)
    .join(vfa1j_ch1)
    .join(vfa2j_ch1)
    .set{vfat1_for_fitting_with_b1}
}else{
vfa1_post_ch1
    .join(vfa_ch1)
    .join(b1_filtered_wo_mask)
    .join(vfa1j_ch1)
    .join(vfa2j_ch1)
    .set{vfat1_for_fitting_with_b1}

}

vfat1_for_fitting_with_b1.into{vfat1_with_b1_bet;vfat1_with_b1}

/*Merge vfa1_post with vfat1_from_alignment only.
WITHOUT B1 MAP
*/
vfa1_post_ch2
    .join(vfa_ch2)
    .join(vfa1j_ch2)
    .join(vfa2j_ch2)
    .set{vfat1_for_fitting_without_b1}

vfat1_for_fitting_without_b1.into{vfat1_without_b1_bet;vfat1_without_b1}

/* We need to join these channels to avoid problems.
WITH B1 MAP
*/
vfat1_with_b1_bet
    .join(mask_from_bet_ch1)
    .set{vfat1_with_b1_bet_merged}

/* Depeding on the nextflow.config 
settings for use_b1cor and use_bet, one of th
following 4 processes will be executed. 
*/

process Fit_VFAT1_With_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true && params.use_bet == true

    input:
        tuple val(sid), file(vfa1), file(vfa2_reg),\
        file(b1map), file(vfa1j), file(vfa2j), file(mask) from vfat1_with_b1_bet_merged
        
    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json"  
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = {'$vfa1, '$vfa2_reg'}; requiredArgs_jsn = {'$vfa1j', '$vfa2j'}; vfa_t1_wrapper2(requiredArgs_nii,requiredArgs_jsn,'mask','$mask','b1map','$b1map','b1factor',$params.b1cor_factor,'qmrlab_path','$params.qmrlab_path', 'sid','${sid}'); exit();"
        """
}

process Fit_VFAT1_With_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_b1cor == true && params.use_bet == false

    input:
        tuple val(sid), file(vfa1), file(vfa2_reg),\
        file(b1map), file(vfa1j), file(vfa2j) from vfat1_with_b1

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz" 
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json" 
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = {'$vfa1, '$vfa2_reg'}; requiredArgs_jsn = {'$vfa1j', '$vfa2j'}; vfa_t1_wrapper2(requiredArgs_nii,requiredArgs_jsn,'b1map','$b1map','b1factor',$params.b1cor_factor,'qmrlab_path','$params.qmrlab_path', 'sid','${sid}'); exit();"
        """
               
}


/* We need to join these channels to avoid problems.*/
vfat1_without_b1_bet
    .join(mask_from_bet_ch2)
    .set{vfat1_without_b1_bet_merged}

process Fit_VFAT1_Without_B1map_With_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == false && params.use_bet==true

    input:
        tuple val(sid), file(vfa1), file(vfa2_reg),\
        file(vfa1j), file(vfa2j), file(mask) from vfat1_without_b1_bet_merged

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz" 
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json" 
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = {'$vfa1, '$vfa2_reg'}; requiredArgs_jsn = {'$vfa1j', '$vfa2j'}; vfa_t1_wrapper2(requiredArgs_nii, requiredArgs_jsn, 'mask','$mask','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"
        """
}

process Fit_VFAT1_Without_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == false && params.use_bet==false

    input:
        tuple val(sid), file(vfa1), file(vfa2_reg),\
        file(vfa1j), file(vfa2j) from vfat1_without_b1

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json"  
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = {'$vfa1, '$vfa2_reg'}; requiredArgs_jsn = {'$vfa1j', '$vfa2j'}; vfa_t1_wrapper2(requiredArgs_nii,requiredArgs_jsn,'qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}
