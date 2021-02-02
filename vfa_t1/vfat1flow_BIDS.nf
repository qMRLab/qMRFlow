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
    in_dataNii = Channel
        .fromFilePairs("$root/**/anat/sub-*_flip-[0-9]_VFA.nii.gz", maxDepth: 4, size: -1, flat: false)
	.transpose()
	.view{"- $it"}

    in_dataNii.into{dataNii_ch1; dataNii_ch2; in_ch3; in_ch4}

    vfa1 = dataNii_ch1
	.first()

    vfa1.into{vfa1_bet; vfa1_alignment}

    vfa1_to_ants = vfa1_alignment
	.flatten()
	.last()
	.view{"- $it"}

    in_dataJSON = Channel
        .fromFilePairs("$root/**/anat/sub-*_flip-[0-9]_VFA.json", maxDepth: 4, size: -1, flat: false)

    /* ==== BIDS: B1 map ==== */             
    /* Look for B1map in fmap folder */
    b1_data = Channel
        .fromFilePairs("$root/**/fmap/sub-*_TB1map.nii.gz", maxDepth:4, size:1, flat:true)
    (b1map) = b1_data       
        .map{sid, B1plusmap -> [tuple(sid, B1plusmap)]}     
        .separate(1)
	      

}   
else{
    error "ERROR: Argument (--root) must be passed. See USAGE."
}

log.info "qMRflow: VFAT1 pipeline"
log.info "======================="
log.info ""
log.info "###   ## ########   ###   ########  ##"
log.info "###   ## ##        ## ##     ##    ###"
log.info "###   ## ##       ##   ##    ##     ##"
log.info "###   ## ####    ##     ##   ##     ##"
log.info "# ## ##  ##      #########   ##     ##"
log.info "#  ###   ##      ##     ##   ##     ##"
log.info "#   #    ##      ##     ##   ##     ##"
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

process Extract_Brain{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    when:
        params.use_bet == true

    input:
        tuple val(sid), file(vfa1) from vfa1_bet

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

/*
dataNii_ch2
    .combine(vfa1_to_ants)
    .view{"- $it"}
*/

process Align_Input_Volumes {
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'

    input:
        tuple val(sid), file(moving) from dataNii_ch2
	each fixed from vfa1_to_ants

    output:
        tuple val(sid), "${moving.simpleName}_aligned.nii.gz"\
        into vfa_aligned
        file "${moving.simpleName}_aligned.nii.gz"
        file "${moving.simpleName}_to_vfa1_displacement.*.mat"

    script:
        """
        antsRegistration -d $params.ants_dim \
                            --float 0 \
                            -o [${moving.simpleName}_to_vfa1_displacement.mat,${moving.simpleName}_aligned.nii.gz] \
                            --transform $params.ants_transform \
                            --metric $params.ants_metric[$fixed,$moving,$params.ants_metric_weight, $params.ants_metric_bins,$params.ants_metric_sampling,$params.ants_metric_samplingprct] \
                            --convergence $params.ants_convergence \
                            --shrink-factors $params.ants_shrink \
                            --smoothing-sigmas $params.ants_smoothing
        """
}

test = vfa_aligned.groupTuple()
	.view{"- $it"}

test2 = in_ch4.groupTuple()
	.view{"- $it"}

process Fit_VFAT1_With_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == true && params.use_bet==false

    input:
        tuple val(sid), file(data1) from test
	tuple val(sid), file(data2) from in_dataJSON
	tuple val(sid), file(b1map) from b1map

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json"  
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = split('$data1'); requiredArgs_jsn = split('$data2'); vfa_t1_wrapper2(requiredArgs_nii',requiredArgs_jsn','b1map','$b1map','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}

/*
process Fit_VFAT1_Without_B1map_Without_Bet{
    tag "${sid}"
    publishDir "$root/derivatives/qMRLab/${sid}", mode: 'copy'
    
    when:
        params.use_b1cor == false && params.use_bet==false

    input:
        tuple val(sid), file(data1) from test2
	tuple val(sid), file(data2) from in_dataJSON

    output:
        file "${sid}_T1map.nii.gz" 
        file "${sid}_M0map.nii.gz"
        file "${sid}_T1map.json" 
        file "${sid}_M0map.json"  
        file "${sid}_vfa_t1.qmrlab.mat"

    script: 
        """
            cp /usr/local/qMRLab/qMRWrappers/vfa_t1/vfa_t1_wrapper2.m vfa_t1_wrapper2.m

            $params.runcmd "requiredArgs_nii = split('$data1'); requiredArgs_jsn = split('$data2'); vfa_t1_wrapper2(requiredArgs_nii',requiredArgs_jsn','qmrlab_path','$params.qmrlab_path', 'sid','${sid}', 'containerType','$workflow.containerEngine', 'containerTag','$params.containerTag', 'description','$params.description', 'datasetDOI','$params.datasetDOI', 'datasetURL','$params.datasetURL', 'datasetVersion','$params.datasetVersion'); exit();"

	    mv dataset_description.json $root/derivatives/qMRLab/dataset_description.json
        """
}
*/



















