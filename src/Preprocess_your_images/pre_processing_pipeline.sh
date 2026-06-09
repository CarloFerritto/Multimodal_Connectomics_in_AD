#!/bin/bash
Machine="cluster"
if [ "$Machine" == "cluster" ]
then
    # Import Anima binaries and scripts
    export PATH=$PATH:/opt/tools/anima/Anima-Binaries-4.2/
    export PATH=$PATH:/opt/tools/anima/Anima-Scripts-Public/
    Anima_dir="/opt/tools/anima/"
    # Import ANTs binaries
    export ANTSPATH=/opt/tools/ants/bin/
    export PATH=${ANTSPATH}:$PATH
    # Import FSL binaries
    FSLDIR=/opt/tools/fsl
    . ${FSLDIR}/etc/fslconf/fsl.sh
    PATH=${FSLDIR}/bin:${PATH}
    export FSLDIR PATH
    ## Activate conda environment
    module load conda
    conda activate project_env
    Dataset="sample_dataset"
    Parcellation_dir="/workspace/project/resources/parcellation/"
    Template_dir="/workspace/project/resources/template/"
    DATASET_dir="/workspace/project/data/${Dataset}/Nifti/"         
fi

Patients=('sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005')


for Patient in "${Patients[@]}" 
do
    cd ${DATASET_dir}


    Cortical_parcellation="Schaefer2018"
    Sub_Cortical_parcellation="Tian_Subcortex_S1"
    Template_Space="MNI152NLin2009cAsym"
    T1_Resolution="1mm"
    N_Cortical_Parcels="400"

    T1_Template_dir=${Template_dir}${Template_Space}/
    T1_Template=${T1_Template_dir}${Template_Space}_T1_${T1_Resolution}.nii.gz
    T1_Template_brain=${T1_Template_dir}${Template_Space}_T1_${T1_Resolution}_brain.nii.gz
    Template_brain_mask=${T1_Template_dir}${Template_Space}_mask.nii
    T1_subject_brain=${Patient}_T1w_masked.nii.gz 
    T1_subject=${Patient}_T1w_reorient.nii.gz
    T1_subject_brain_mask=${Patient}_T1w_brainMask.nii.gz  
    T1_wmseg=${Patient}_T1w_masked_wmseg.nii.gz

    echo "========================================"
    echo "Processing patient:" ${Patient}

    # Preprocessing variables 
    Manufacturer=$(jq  '.Manufacturer' ${Patient}/anat/${Patient}_T1w.json)
    if [[ "$Manufacturer" == *"Philips"* ]]
        then
        dwelltime_dwi=$(jq  '.EstimatedEffectiveEchoSpacing' ${Patient}/dwi/${Patient}_dwi.json)
        readout_time_dwi=$(jq  '.EstimatedTotalReadoutTime' ${Patient}/dwi/${Patient}_dwi.json)
        phase_encoding_dwi=$(jq  '.PhaseEncodingAxis' ${Patient}/dwi/${Patient}_dwi.json)
        dwelltime_fmri=$(jq  '.EstimatedEffectiveEchoSpacing' ${Patient}/func/${Patient}_task-rest_bold.json)
        readout_time_fmri=$(jq  '.EstimatedTotalReadoutTime' ${Patient}/func/${Patient}_task-rest_bold.json)
        TR_fmri=$(jq  '.RepetitionTime' ${Patient}/func/${Patient}_task-rest_bold.json)
        phase_encoding_fmri=$(jq  '.PhaseEncodingAxis' ${Patient}/func/${Patient}_task-rest_bold.json)
        echo1_fmap=$(jq  '.EchoTime' ${Patient}/fmap/${Patient}_echo-1_part-mag.json)
        echo2_fmap=$(jq  '.EchoTime' ${Patient}/fmap/${Patient}_echo-2_part-mag.json)
    else
        dwelltime_dwi=$(jq  '.EffectiveEchoSpacing' ${Patient}/dwi/${Patient}_dwi.json)
        readout_time_dwi=$(jq  '.TotalReadoutTime' ${Patient}/dwi/${Patient}_dwi.json)
        phase_encoding_dwi=$(jq  '.PhaseEncodingDirection' ${Patient}/dwi/${Patient}_dwi.json)
        dwelltime_fmri=$(jq  '.EffectiveEchoSpacing' ${Patient}/func/${Patient}_task-rest_bold.json)
        readout_time_fmri=$(jq  '.TotalReadoutTime' ${Patient}/func/${Patient}_task-rest_bold.json)
        TR_fmri=$(jq  '.RepetitionTime' ${Patient}/func/${Patient}_task-rest_bold.json)
        slice_timing_fmri=$(jq  '.SliceTiming' ${Patient}/func/${Patient}_task-rest_bold.json)
        phase_encoding_fmri=$(jq  '.PhaseEncodingDirection' ${Patient}/func/${Patient}_task-rest_bold.json)
        echo1_fmap=$(jq  '.EchoTime' ${Patient}/fmap/${Patient}_echo-1_part-mag.json)
        echo2_fmap=$(jq  '.EchoTime' ${Patient}/fmap/${Patient}_echo-2_part-mag.json)
    fi
    python3 metadata_handler.py -dir $DATASET_dir -p $Patient -m $Manufacturer -dd $dwelltime_dwi -rotd $readout_time_dwi -ped $phase_encoding_dwi -df $dwelltime_fmri -rotf $readout_time_fmri -trf $TR_fmri -pef $phase_encoding_fmri 
    
    cd ${DATASET_dir}${Patient}/anat/
    
    echo "========================================"
    echo "Processing T1"
    # Bias Field Correction [ANTS: N4BiasFieldCorrection]
    N4BiasFieldCorrection -i ${Patient}_T1w.nii.gz -o ${Patient}_T1w_bias.nii.gz
    
    # Reorientation (rotations of 0, 90, 180, 270° only)
    echo " "
    echo " => Reorientation [fslreorient2std]"
    fslreorient2std ${Patient}_T1w_bias.nii.gz ${Patient}_T1w_reorient.nii.gz

    # Brain extraction and masking
    echo " "
    echo " => Brain extraction and masking [animaAtlasBasedBrainExtraction.py]"
    python3 ${Anima_dir}/Anima-Scripts-Public/brain_extraction/animaAtlasBasedBrainExtraction.py -i ${Patient}_T1w_reorient.nii.gz
    animaConvertImage -i ${Patient}_T1w_reorient_brainMask.nrrd -o ${Patient}_T1w_brainMask.nii.gz
    animaConvertImage -i ${Patient}_T1w_reorient_masked.nrrd -o ${Patient}_T1w_masked.nii.gz

    # Brain segmentation (WM/GM/CSF) using FAST on the masked brain
    echo " "
    echo " => Brain segmentation (WM/GM/CSF) in subject's space [FSL: fast]"
    fast -o ${Patient}_T1w_masked ${Patient}_T1w_masked.nii.gz
    fslmaths ${Patient}_T1w_masked_pve_2 -thr 0.5 -bin ${Patient}_T1w_masked_wmseg
    
    # Set parcellation name
    
    if [[ "$Cortical_parcellation" == "" ]]  &&  [[ "$Sub_Cortical_parcellation" == "" ]]
    then
        exit "Please provide either a cortical or subcortical parcellation!"
    elif [[ "$Sub_Cortical_parcellation" == "" ]]
    then
        echo "Using only $Cortical_parcellation cortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Cortex/${Cortical_parcellation}_${N_Cortical_Parcels}Parcels_7Networks_order_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks
    elif [[ "Cortical_parcellation" == "" ]]
    then
        echo "Using only $Sub_Cortical_parcellation subcortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Subcortex/${Sub_Cortical_parcellation}_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Sub_Cortical_parcellation}
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Sub_Cortical_parcellation}
    else    
        echo "Using $Cortical_parcellation cortical parcellation and $Sub_Cortical_parcellation subcortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Cortex-Subcortex/${Cortical_parcellation}_${N_Cortical_Parcels}Parcels_7Networks_order_${Sub_Cortical_parcellation}_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}
    fi 

    # Use ANTs for nonlinear registration to template space
    echo " =>"
    echo "Nonlinear registration to template space ($Template_Space) [ANTS: antsRegistration (Rigid+Affine+NonLinear)]"
    antsRegistration -d 3 --float 0 -o [ ${Patient}_T1to${Template_Space}_ , ${Patient}_T1to${Template_Space}.nii.gz] -n Linear -w [ 0.005 , 0.995] -u 0 -r [ $T1_Template_brain, $T1_subject_brain, 1] -t Rigid[0.1] -m MI[$T1_Template_brain, $T1_subject_brain, 1, 32, Regular, 0.25] -c [ 1000x500x250x100, 1e-7, 10]  -f 8x4x2x1 -s 3x2x1x0vox -t Affine[0.1] -m MI[ $T1_Template_brain, $T1_subject_brain, 1, 32, Regular, 0.25] -c [ 1000x500x250x100, 1e-7,10] -f 8x4x2x1 -s 3x2x1x0vox -t SyN[ 0.1, 3, 0] -m CC[ $T1_Template_brain, $T1_subject_brain, 1, 4] -c [ 200x200x200x200, 1e-7, 10] -f 8x4x2x1 -s 3x2x1x0vox
    
    echo " =>"
    echo "Registering parcellation to subject space [ANTS: antsApplyTransforms]"
    antsApplyTransforms -d 3 --float 0 -i $Parcellation  -r $T1_subject_brain -n NearestNeighbor -t [ ${Patient}_T1to${Template_Space}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space}_1InverseWarp.nii.gz  -o ${Parcellation_T1_space}.nii.gz -v 1
    
    echo " =>"
    echo "Registering template mask to subject space [ANTS: antsApplyTransforms]"
    antsApplyTransforms -d 3 --float 0 -i $Template_brain_mask  -r $T1_subject_brain -n NearestNeighbor -t [ ${Patient}_T1to${Template_Space}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space}_1InverseWarp.nii.gz  -o ${Patient}_T1w_brainMask_tight.nii.gz -v 1
    fslmaths $T1_subject_brain_mask -mul ${Patient}_T1w_brainMask_tight.nii.gz ${Patient}_T1w_brainMask_tight.nii.gz

    
    echo "========================================"
    echo "Processing diffusion MRI"

    cd ${DATASET_dir}${Patient}/dwi/ 
    # Reorient DWI volumes
    dwi_file=${Patient}_dwi
    fslreorient2std ${dwi_file}.nii.gz  ${dwi_file}_reorient.nii.gz 
    # Estimate masks from T1
    echo " =>"
    echo "DWI mask estimation"
    dwi_mask=${Patient}_T1w_brainMask_resampled_to_dwi_space
    B0_file=${Patient}_dwi_B0_raw
    B0_mean=${Patient}_dwi_B0_mean_raw
    
    echo " =>"
    echo "Selecting B0 volumes"
    dwiextract -force -bzero -fslgrad ${dwi_file}_real.bvec ${dwi_file}.bval ${dwi_file}_reorient.nii.gz ${B0_file}.nii.gz
    mrmath ${B0_file}.nii.gz mean ${B0_mean}.nii.gz -axis 3 -force
    fslroi ${B0_file}.nii.gz ${B0_file}_first.nii.gz 0 1
    antsRegistration -d 3 --float 0 -o [ ${Patient}_T1toDWI_ , ${Patient}_T1toDWI.nii.gz] -n Linear -w [ 0.005 , 0.995] -u 0 -r [ ${B0_file}_first.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked.nii.gz, 1] -t Rigid[0.1] -m MI[${B0_mean}.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked.nii.gz, 1, 32, Regular, 0.25] -c [ 1000x1000x1000x1000, 1e-7, 10]  -f 8x4x2x1 -s 3x2x1x0vox 
    antsApplyTransforms -d 3 --float 0 -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask.nii.gz -r ${B0_file}_first.nii.gz -n NearestNeighbor -t ${Patient}_T1toDWI_0GenericAffine.mat -o ${dwi_mask}.nii.gz

    # Extract the number of volumes
    nvol=$(fslnvols ${dwi_file}_reorient.nii.gz)
    # Create a dummy index file (direction is always PA) for eddy
    index=""
    for ((j=1; j<=nvol; j+=1)); do index="$index 1"; done
    dwi_txt2=${dwi_file}_index.txt
    echo $index > $dwi_txt2

    dwi_eddy=${dwi_file}_eddy
    echo " =>"
    echo "Using eddy without fieldmap"
    eddy --imain=${dwi_file}_reorient.nii.gz --mask=${dwi_mask}.nii.gz --acqp=${DATASET_dir}${Patient}/config_dwi.txt --index=$dwi_txt2 --bvecs=${dwi_file}_real.bvec --bvals=${dwi_file}.bval --repol --out=$dwi_eddy --data_is_shelled

    if [[ "$Manufacturer" == *"Siemens"* ]]  || [[ "$Manufacturer" == *"Philips"* ]]
    then
        if [[ "$Manufacturer" == *"Siemens"* ]]
        then
            ### FIELDMAP PREPARATION SIEMENS
            echo " =>"
            echo "Preparing fieldmap"
            cd ${DATASET_dir}${Patient}/fmap/
            field_mag=${Patient}_echo-2_part-mag
            field_phase=${Patient}_phasediff
            field_mag_masked=${Patient}_magnitude_masked
            field_mag_mask=${Patient}_magnitude_brainMask  
            field_map_rad=${Patient}_fieldmap_rad
            field_map_hz=${Patient}_fieldmap
            echotime_diff=$(bc <<< "scale=10; $echo2_fmap-$echo1_fmap")
            echotime_diff=$(bc <<< "scale=10; $echotime_diff*1000")
            fslreorient2std $field_mag ${field_mag}_reorient
            fslreorient2std $field_phase ${field_phase}_reorient

            antsRegistration -d 3 --float 0 -o [ ${Patient}_T1tofieldmap_ , ${Patient}_T1tofieldmap.nii.gz] -n Linear -w [ 0.005 , 0.995] -u 0 -r [ ${field_mag}_reorient.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_reorient.nii.gz, 1] -t Rigid[0.1] -m MI[${field_mag}_reorient.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_reorient.nii.gz, 1, 32, Regular, 0.25] -c [ 1000x1000x1000x1000, 1e-7, 10]  -f 8x4x2x1 -s 3x2x1x0vox 
            antsApplyTransforms -d 3 --float 0 -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_tight.nii.gz -r ${field_mag}_reorient.nii.gz -n Linear -t ${Patient}_T1tofieldmap_0GenericAffine.mat -o ${field_mag_mask}_linear.nii.gz
            fslmaths ${field_mag_mask}_linear.nii.gz -thr 0.999 -bin ${field_mag_mask}.nii.gz
            fslmaths ${field_mag}_reorient -mul ${field_mag_mask}.nii.gz $field_mag_masked

            fsl_prepare_fieldmap SIEMENS ${field_phase}_reorient ${field_mag_masked} $field_map_rad $echotime_diff
            fugue --loadfmap=$field_map_rad --despike -m --savefmap=$field_map_rad 

        elif [[ "$Manufacturer" == *"Philips"* ]]
        then
            ### FIELDMAP PREPARATION PHILIPS
            echo " =>"
            echo "Preparing fieldmap"
            cd ${DATASET_dir}${Patient}/fmap/
            field_mag=${Patient}_echo-2_part-mag
            field_phase_1=${Patient}_echo-1_part-phase
            field_phase_2=${Patient}_echo-2_part-phase
            field_mag_masked=${Patient}_magnitude_masked
            field_mag_mask=${Patient}_magnitude_brainMask  
            field_map_rad=${Patient}_fieldmap_rad
            echotime_diff=$(bc <<< "scale=10; $echo2_fmap-$echo1_fmap")
            echotime_diff=$(bc <<< "scale=10; $echotime_diff*1000")
            fslreorient2std $field_mag ${field_mag}_reorient
            fslreorient2std $field_phase_1 ${field_phase_1}_reorient
            fslreorient2std $field_phase_2 ${field_phase_2}_reorient

            fslmaths ${field_phase_2}_reorient -add 3.14159 ${field_phase_2}_reorient_rad -odt float
            fslmaths ${field_phase_1}_reorient -add 3.14159 ${field_phase_1}_reorient_rad -odt float
            antsRegistration -d 3 --float 0 -o [ ${Patient}_T1tofieldmap_ , ${Patient}_T1tofieldmap.nii.gz] -n Linear -w [ 0.005 , 0.995] -u 0 -r [ ${field_mag}_reorient.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_reorient.nii.gz, 1] -t Rigid[0.1] -m MI[${field_mag}_reorient.nii.gz, ${DATASET_dir}${Patient}/anat/${Patient}_T1w_reorient.nii.gz, 1, 32, Regular, 0.25] -c [ 1000x1000x1000x1000, 1e-7, 10]  -f 8x4x2x1 -s 3x2x1x0vox 
            antsApplyTransforms -d 3 --float 0 -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_tight.nii.gz -r ${field_mag}_reorient.nii.gz -n Linear -t ${Patient}_T1tofieldmap_0GenericAffine.mat -o ${field_mag_mask}_linear.nii.gz
            fslmaths ${field_mag_mask}_linear.nii.gz -thr 0.999 -bin ${field_mag_mask}.nii.gz
            fslmaths ${field_mag}_reorient -mul ${field_mag_mask}.nii.gz $field_mag_masked

            prelude -a ${field_mag}_reorient -p ${field_phase_1}_reorient_rad.nii.gz -o ${field_phase_1}_reorient_rad_unwarp.nii.gz
            prelude -a ${field_mag}_reorient -p ${field_phase_2}_reorient_rad.nii.gz -o ${field_phase_2}_reorient_rad_unwarp.nii.gz
            fslmaths ${field_phase_2}_reorient_rad_unwarp.nii.gz -sub ${field_phase_1}_reorient_rad_unwarp.nii.gz -mul 1000 -div $echotime_diff $field_map_rad -odt float
            fslmaths $field_map_rad -mul ${field_mag_mask}.nii.gz $field_map_rad
            fugue --loadfmap=$field_map_rad --despike -m --savefmap=$field_map_rad
        fi
        ### Correct for distortion
        cd ${DATASET_dir}${Patient}/dwi/
        dwiextract -force -bzero -fslgrad ${dwi_file}_real.bvec ${dwi_file}.bval ${dwi_eddy}.nii.gz ${B0_file}_eddy.nii.gz
        fslroi ${B0_file}_eddy.nii.gz ${B0_file}_eddy_first.nii.gz 0 1
        #BBR Registration
        echo " =>"
        echo "BBR registration to T1"
        if [[ "$phase_encoding_dwi" == *"j"* ]] && [[ "$phase_encoding_dwi" != *"-"* ]]
        then
            echo "direction y"
            epi_reg --epi=${B0_file}_eddy_first.nii.gz --t1=${DATASET_dir}${Patient}/anat/${T1_subject} --t1brain=${DATASET_dir}${Patient}/anat/${T1_subject_brain} --wmseg=${DATASET_dir}${Patient}/anat/${T1_wmseg} --out=${Patient}_DWI2T1 --fmap=${DATASET_dir}${Patient}/fmap/$field_map_rad  --fmapmag=${DATASET_dir}${Patient}/fmap/${field_mag}_reorient --fmapmagbrain=${DATASET_dir}${Patient}/fmap/$field_mag_masked --echospacing=$dwelltime_dwi --pedir=y
            # Apply FUGUE
            echo " =>"
            echo "Apply FUGUE"
            fugue -i ${dwi_eddy}.nii.gz --dwell=$dwelltime_dwi --loadfmap=${Patient}_DWI2T1_fieldmaprads2epi.nii.gz  --unwarpdir=y -u ${dwi_eddy}_undistorted.nii.gz
        fi
        if [[ "$phase_encoding_dwi" == *"j-"* ]]
        then
            echo "direction y-"
            epi_reg --epi=${B0_file}_eddy_first.nii.gz --t1=${DATASET_dir}${Patient}/anat/${T1_subject} --t1brain=${DATASET_dir}${Patient}/anat/${T1_subject_brain} --wmseg=${DATASET_dir}${Patient}/anat/${T1_wmseg} --out=${Patient}_DWI2T1 --fmap=${DATASET_dir}${Patient}/fmap/$field_map_rad  --fmapmag=${DATASET_dir}${Patient}/fmap/${field_mag}_reorient --fmapmagbrain=${DATASET_dir}${Patient}/fmap/$field_mag_masked --echospacing=$dwelltime_dwi --pedir=y-
            # Apply FUGUE
            echo " =>"
            echo "Apply FUGUE"
            fugue -i ${dwi_eddy}.nii.gz --dwell=$dwelltime_dwi --loadfmap=${Patient}_DWI2T1_fieldmaprads2epi.nii.gz  --unwarpdir=y- -u ${dwi_eddy}_undistorted.nii.gz
        fi
        if [[ "$phase_encoding_dwi" == "i" ]] && [[ "$phase_encoding_dwi" != *"-"* ]]
        then
            echo "direction x"
        fi
        if [[ "$phase_encoding_dwi" == "i-" ]]
        then
            echo "direction x-"
        fi
        if [[ "$phase_encoding_dwi" == "k" ]] && [[ "$phase_encoding_dwi" != *"-"* ]]
        then
            echo "direction z"
        fi
        if [[ "$phase_encoding_dwi" == "k-" ]]
        then
            echo "direction z-"
        fi

    elif [[ "$Manufacturer" == *"GE"* ]]
    then
        ### Correct distortion with Anima
        dwiextract -force -bzero -fslgrad ${dwi_file}_real.bvec ${dwi_file}.bval ${dwi_eddy}.nii.gz ${B0_file}_eddy.nii.gz
        fslroi ${B0_file}_eddy.nii.gz ${B0_file}_eddy_first.nii.gz 0 1
        animaConvertImage -i ${B0_file}_eddy_first.nii.gz  -o ${B0_file}_eddy_first.nrrd
        animaConvertImage -i ${Patient}_T1toDWI.nii.gz  -o ${Patient}_T1toDWI.nrrd
        animaConvertImage -i ${dwi_eddy}.nii.gz -o ${dwi_eddy}.nrrd

        
        if [[ "$phase_encoding_dwi" == *"i"* ]]
        then
            echo "direction x"
            animaDenseSVFBMRegistration -r ${Patient}_T1toDWI.nrrd -m ${B0_file}_eddy_first.nrrd -o ${B0_file}_eddy_first_corrected.nrrd -d 0 -O ${B0_file}_correction_tr.nrrd -t 3 --sym-reg 2
        elif [[ "$phase_encoding_dwi" == *"j"* ]]
        then
            echo "direction y"
            animaDenseSVFBMRegistration -r ${Patient}_T1toDWI.nrrd -m ${B0_file}_eddy_first.nrrd -o ${B0_file}_eddy_first_corrected.nrrd -d 1 -O ${B0_file}_correction_tr.nrrd -t 3 --sym-reg 2
        elif [[ "$phase_encoding_dwi" == *"k"* ]]
        then
            echo "direction z"
            animaDenseSVFBMRegistration -r ${Patient}_T1toDWI.nrrd -m ${B0_file}_eddy_first.nrrd -o ${B0_file}_eddy_first_corrected.nrrd -d 2 -O ${B0_file}_correction_tr.nrrd -t 3 --sym-reg 2
        fi
        animaApplyDistortionCorrection -f ${dwi_eddy}.nrrd -t ${B0_file}_correction_tr.nrrd -o ${dwi_eddy}_undistorted.nii.gz

        rm ${dwi_eddy}.nrrd

    else  
        exit  "$Manufacturer is not a supported manufacturer. Please check the input metadata"
    fi 

    echo " "
    echo " => Bias field correction on DWI"
    dwibiascorrect ants ${dwi_eddy}_undistorted.nii.gz ${dwi_eddy}_undistorted.nii.gz -fslgrad ${dwi_file}_real.bvec ${dwi_file}.bval  -force
    
    dwiextract -force -bzero -fslgrad ${dwi_file}_real.bvec ${dwi_file}.bval ${dwi_eddy}_undistorted.nii.gz ${B0_file}_eddy_undistorted.nii.gz
    fslroi ${B0_file}_eddy_undistorted.nii.gz ${B0_file}_eddy_undistorted_first.nii.gz 0 1
    echo "BBR registration to T1 after fugue"
    flirt -ref ${DATASET_dir}${Patient}/anat/${T1_subject} -in ${B0_file}_eddy_undistorted_first.nii.gz  -dof 6 -omat ${Patient}_DWI2T1_init.mat
    flirt -ref ${DATASET_dir}${Patient}/anat/${T1_subject} -in ${B0_file}_eddy_undistorted_first.nii.gz  -dof 6 -cost bbr -wmseg ${DATASET_dir}${Patient}/anat/${T1_wmseg} -init ${Patient}_DWI2T1_init.mat -omat ${Patient}_DWI2T1.mat -out ${Patient}_DWI2T1.nii.gz -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
    convert_xfm -omat ${Patient}_DWI2T1_inv.mat -inverse ${Patient}_DWI2T1.mat
    reference_image=${B0_file}_eddy_undistorted_first.nii.gz 

    dwi_fin=${dwi_file}_preprocessed
    cp ${dwi_file}_eddy.eddy_rotated_bvecs ${dwi_fin}.bvec
    cp ${dwi_file}.bval ${dwi_fin}.bval 

    # Reorient T1 in diffusion space
    echo " "
    echo " => Reorienting T1 in diffusion space based on BBR registration and without resampling"
    flirt  -in ${DATASET_dir}${Patient}/anat/${T1_subject_brain_mask} -ref $reference_image -interp nearestneighbour -applyxfm -init ${Patient}_DWI2T1_inv.mat -out ${Patient}_dwi_brainMask.nii.gz
    transformconvert ${Patient}_DWI2T1.mat $reference_image  ${DATASET_dir}${Patient}/anat/${T1_subject_brain} flirt_import ${Patient}_DWI2T1_mrtrix.txt -force
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse ${DATASET_dir}${Patient}/anat/${T1_subject_brain} ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg.nii.gz -force
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse -interp nearest ${Parcellation_T1_space}.nii.gz ${Parcellation_T1_space}_coreg.nii.gz -force

    
    # Re-orient images to be axial first
    echo " "
    echo " => Reorienting images"
    animaConvertImage -i ${dwi_eddy}_undistorted.nii.gz -o ${dwi_fin}_or.nrrd -R AXIAL
    animaConvertImage -i ${Patient}_dwi_brainMask.nii.gz -o ${Patient}_dwi_brainMask_axial.nrrd -R AXIAL
    animaConvertImage -i ${Parcellation_T1_space}_coreg.nii.gz -o ${Parcellation_T1_space}_coreg_axial.nrrd -R AXIAL
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg.nii.gz  -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial.nrrd -R AXIAL
    
    # Perform denoising
    echo " "
    echo " => Performing denoising"
    animaNLMeansTemporal -i ${dwi_fin}_or.nrrd -b 0.5 -n 3 -o ${dwi_fin}_or_nlm.nrrd
    animaConvertImage -i ${dwi_fin}_or_nlm.nrrd  -o ${dwi_fin}_or_nlm.nii.gz
    animaConvertImage -i ${Patient}_dwi_brainMask_axial.nrrd  -o ${Patient}_dwi_brainMask_axial.nii.gz
    animaConvertImage -i ${Parcellation_T1_space}_coreg_axial.nrrd -o ${Parcellation_T1_space}_coreg_axial.nii.gz
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial.nrrd -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial.nii.gz
    
    # Brain masking
    echo " "
    echo " => Brain masking image"

    B0_file=${Patient}_dwi_B0_skull
    B0_mean_file=${Patient}_dwi_B0_skull_mean
    dwiextract -force -bzero -fslgrad ${dwi_fin}.bvec ${dwi_fin}.bval ${dwi_fin}_or_nlm.nii.gz ${B0_file}.nii.gz
    mrmath ${B0_file}.nii.gz mean ${B0_mean_file}.nii.gz -axis 3 -force
    fslmaths ${dwi_fin}_or_nlm.nii.gz -mul ${Patient}_dwi_brainMask_axial.nii.gz ${dwi_fin}.nii.gz

    B0_file=${Patient}_dwi_B0
    B0_mean_file=${Patient}_dwi_B0_mean
    dwiextract -force -bzero -fslgrad ${dwi_fin}.bvec ${dwi_fin}.bval ${dwi_fin}.nii.gz ${B0_file}.nii.gz
    mrmath ${B0_file}.nii.gz mean ${B0_mean_file}.nii.gz -axis 3 -force

    # Estimate tensors
    Tensors_file=${dwi_fin}_tensors
    animaConvertImage -i ${dwi_fin}.nii.gz  -o ${dwi_fin}.nrrd
    animaDTIEstimator -i ${dwi_fin}.nrrd -o ${dwi_fin}_Tensors.nrrd -g ${dwi_fin}.bvec -b ${dwi_fin}.bval
    rm ${dwi_fin}.nrrd
    rm ${dwi_fin}_or_nlm.nrrd
    rm ${dwi_fin}_or.nrrd

    echo "========================================"
    echo "Processing fMRI"
    cd ${DATASET_dir}${Patient}/func/
    fMRI_file=${Patient}_task-rest_bold
    Parcellation_fMRI_space=${DATASET_dir}${Patient}/func/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}
    
    fslreorient2std ${fMRI_file}.nii.gz  ${fMRI_file}_reorient.nii.gz
    
    echo " "
    echo " => Motion correction [FSL:mcflirt]"
    mcflirt -in ${fMRI_file}_reorient.nii.gz -out ${fMRI_file}_motion_corrected.nii.gz -dof 6 -mats -plots # rBOLD_reorient.nii.gz

    echo " "
    echo " => Slice timing (please modify the −−odd option if you have different order) [FSL:slicetimer]"
    slicetimer -i ${fMRI_file}_motion_corrected.nii.gz -o ${fMRI_file}_slicetimer.nii.gz -r ${TR_fmri} --tcustom=${DATASET_dir}${Patient}/sliceTimer.txt -v 
    fslroi ${fMRI_file}_slicetimer.nii.gz ${fMRI_file}_slicetimer_first.nii.gz 0 1

    if [[ "$Manufacturer" == *"Siemens"* ]]  || [[ "$Manufacturer" == *"Philips"* ]]
    then

        if [[ "$phase_encoding_fmri" == *"j"* ]] && [[ "$phase_encoding_fmri" != *"-"* ]]
        then
            echo "direction y"
            epi_reg --epi=${fMRI_file}_slicetimer_first.nii.gz --t1=${DATASET_dir}${Patient}/anat/${T1_subject} --t1brain=${DATASET_dir}${Patient}/anat/${T1_subject_brain} --wmseg=${DATASET_dir}${Patient}/anat/${T1_wmseg} --out=${Patient}_FMRI2T1 --fmap=${DATASET_dir}${Patient}/fmap/$field_map_rad  --fmapmag=${DATASET_dir}${Patient}/fmap/${field_mag}_reorient --fmapmagbrain=${DATASET_dir}${Patient}/fmap/$field_mag_masked --echospacing=$dwelltime_fmri --pedir=y
            # Apply FUGUE
            echo " =>"
            echo "Apply FUGUE"
            fugue -i ${fMRI_file}_slicetimer.nii.gz --dwell=$dwelltime_fmri --loadfmap=${Patient}_FMRI2T1_fieldmaprads2epi.nii.gz  --unwarpdir=y -u ${fMRI_file}_disto.nii.gz
        fi
        if [[ "$phase_encoding_fmri" == *"j-"* ]]
        then
            echo "direction y-"
            epi_reg --epi=${fMRI_file}_slicetimer_first.nii.gz --t1=${DATASET_dir}${Patient}/anat/${T1_subject} --t1brain=${DATASET_dir}${Patient}/anat/${T1_subject_brain} --wmseg=${DATASET_dir}${Patient}/anat/${T1_wmseg} --out=${Patient}_FMRI2T1 --fmap=${DATASET_dir}${Patient}/fmap/$field_map_rad  --fmapmag=${DATASET_dir}${Patient}/fmap/${field_mag}_reorient --fmapmagbrain=${DATASET_dir}${Patient}/fmap/$field_mag_masked --echospacing=$dwelltime_fmri --pedir=y
            # Apply FUGUE
            echo " =>"
            echo "Apply FUGUE"
            fugue -i ${fMRI_file}_slicetimer.nii.gz --dwell=$dwelltime_fmri --loadfmap=${Patient}_FMRI2T1_fieldmaprads2epi.nii.gz  --unwarpdir=y- -u ${fMRI_file}_disto.nii.gz
        fi
        if [[ "$phase_encoding_fmri" == "i" ]] && [[ "$phase_encoding_fmri" != *"-"* ]]
        then
            echo "direction x"
        fi
        if [[ "$phase_encoding_fmri" == "i-" ]]
        then
            echo "direction x-"
        fi
        if [[ "$phase_encoding_fmri" == "k" ]] && [[ "$phase_encoding_fmri" != *"-"* ]]
        then
            echo "direction z"
        fi
        if [[ "$phase_encoding_fmri" == "k-" ]]
        then
            echo "direction z-"
        fi

    elif [[ "$Manufacturer" == *"GE"* ]]
    then
        #BBR Registration
        echo "BBR registration to T1"
        epi_reg --epi=${fMRI_file}_slicetimer_first.nii.gz --t1=${DATASET_dir}${Patient}/anat/${T1_subject} --t1brain=${DATASET_dir}${Patient}/anat/${T1_subject_brain} --wmseg=${DATASET_dir}${Patient}/anat/${T1_wmseg} --out=${Patient}_FMRI2T1

        ##Register T1 to fMRI space 
        convert_xfm -omat ${Patient}_T12FMRI.mat -inverse ${Patient}_FMRI2T1.mat
        flirt  -in ${DATASET_dir}${Patient}/anat/${T1_subject_brain} -ref ${fMRI_file}_slicetimer_first.nii.gz -applyxfm -init ${Patient}_T12FMRI.mat -out ${Patient}_T12FMRI.nii.gz
        ### Correct distortion with Anima
        animaConvertImage -i ${fMRI_file}_slicetimer_first.nii.gz  -o ${fMRI_file}_slicetimer_first.nrrd
        animaConvertImage -i ${Patient}_T12FMRI.nii.gz  -o ${Patient}_T12FMRI.nrrd
        animaConvertImage -i ${fMRI_file}_slicetimer.nii.gz -o ${fMRI_file}_slicetimer.nrrd

        if [[ "$phase_encoding_fmri" == *"i"* ]]
        then
            echo "direction x"
            animaDenseSVFBMRegistration -r ${Patient}_T12FMRI.nrrd -m ${fMRI_file}_slicetimer_first.nrrd -o ${fMRI_file}_slicetimer_first_corrected.nrrd -d 0 -O ${fMRI_file}_correction_tr.nrrd -t 3 --sym-reg 2
        elif [[ "$phase_encoding_fmri" == *"j"* ]]
        then
            echo "direction y"
            animaDenseSVFBMRegistration -r ${Patient}_T12FMRI.nrrd -m ${fMRI_file}_slicetimer_first.nrrd -o ${fMRI_file}_slicetimer_first_corrected.nrrd -d 1 -O ${fMRI_file}_correction_tr.nrrd -t 3 --sym-reg 2
        elif [[ "$phase_encoding_fmri" == *"k"* ]]
        then
            echo "direction z"
            animaDenseSVFBMRegistration -r ${Patient}_T12FMRI.nrrd -m ${fMRI_file}_slicetimer_first.nrrd -o ${fMRI_file}_slicetimer_first_corrected.nrrd -d 2 -O ${fMRI_file}_correction_tr.nrrd -t 3 --sym-reg 2
        fi
        animaApplyDistortionCorrection -f ${fMRI_file}_slicetimer.nrrd -t ${fMRI_file}_correction_tr.nrrd -o ${fMRI_file}_disto.nrrd
        animaConvertImage -i ${fMRI_file}_disto.nrrd -o ${fMRI_file}_disto.nii.gz
        rm ${fMRI_file}_disto.nrrd
        rm ${fMRI_file}_slicetimer.nrrd
    else  
        exit  "$Manufacturer is not a supported manufacturer. Please check the input metadata"
    fi 

    
    fslroi ${fMRI_file}_disto.nii.gz ${fMRI_file}_disto_first.nii.gz 0 1

    echo "BBR registration to T1 after fugue"
    flirt -ref ${DATASET_dir}${Patient}/anat/${T1_subject} -in ${fMRI_file}_disto_first.nii.gz -dof 6 -omat ${Patient}_FMRI2T1_init.mat
    flirt -ref ${DATASET_dir}${Patient}/anat/${T1_subject} -in ${fMRI_file}_disto_first.nii.gz -dof 6 -cost bbr -wmseg ${DATASET_dir}${Patient}/anat/${T1_wmseg} -init ${Patient}_FMRI2T1_init.mat -omat ${Patient}_FMRI2T1_fin.mat -out ${Patient}_FMRI2T1_fin.nii.gz -schedule ${FSLDIR}/etc/flirtsch/bbr.sch
    convert_xfm -omat ${Patient}_T12FMRI_fin.mat -inverse ${Patient}_FMRI2T1_fin.mat
    flirt -in ${DATASET_dir}${Patient}/anat/${T1_subject_brain_mask} -ref ${fMRI_file}_disto_first.nii.gz  -interp nearestneighbour -applyxfm -init ${Patient}_T12FMRI_fin.mat -out ${Patient}_fMRI_brainMask.nii.gz
    flirt -in ${Parcellation_T1_space}.nii.gz -ref ${fMRI_file}_disto_first.nii.gz  -interp nearestneighbour -applyxfm -init ${Patient}_T12FMRI_fin.mat -out ${Parcellation_fMRI_space}.nii.gz
    flirt -in ${DATASET_dir}${Patient}/anat/${T1_subject_brain} -ref ${fMRI_file}_disto_first.nii.gz -applyxfm -init ${Patient}_T12FMRI_fin.mat -out ${Patient}_T12FMRI_fin.nii.gz

    echo " Median rBOLD volume"
    fslmaths ${fMRI_file}_disto.nii.gz -Tmedian ${fMRI_file}_disto_med.nii.gz

    echo " => Global intensity normalisation [FSL:fslmaths]"
    fslmaths ${fMRI_file}_disto.nii.gz -ing 10000 ${fMRI_file}_norm.nii.gz -odt float
    cd ${DATASET_dir}${Patient}/anat/

    echo "========================================"
    echo "Tractography"
    T1=${Patient}_T1w_masked_coreg_axial
    fivett_file=${DATASET_dir}${Patient}/anat/${T1}_5tt
    seeding_mask=${DATASET_dir}${Patient}/anat/${T1}_5tt2gmwmi

    echo " "
    echo " => Generating 5tt image from T1 at 1mm resolution "
    5ttgen fsl ${T1}.nii.gz ${fivett_file}.nii.gz -premasked -nocrop -force
    5tt2gmwmi ${fivett_file}.nii.gz ${seeding_mask}.nii.gz -force #new

    echo "========================================"
    echo "Processing DWI for mrtrix space"

    cd ${DATASET_dir}${Patient}/dwi/
    dwi_file=${Patient}_dwi_preprocessed
    
    echo " "
    echo " => Converting bvec file mrtrix space coordinates "
        if [ "$Machine" == "cluster" ]
    then 
        python3 /home/user/shared-storage/share/users/user/flip_bvec.py -i ${dwi_file}.bvec -d 0 -o ${dwi_file}_mrtrix.bvec
    elif [ "$Machine" == "cluster" ]
    then
        python3 /home/usero/Scripts/Preprocessing/diffusion/Tractography/flip_bvec.py -i ${dwi_file}.bvec -d 0 -o ${dwi_file}_mrtrix.bvec
    fi
    

    echo " "
    echo " => Response function "
    dwi2response dhollander ${dwi_file}.nii.gz  ${dwi_file}_out_wm ${dwi_file}_out_gm ${dwi_file}_out_csf -fslgrad ${dwi_file}_mrtrix.bvec ${dwi_file}.bval -force -shells 0,1000 

    echo " "
    echo " => FOD "
    dwi2fod msmt_csd ${dwi_file}.nii.gz ${dwi_file}_out_wm ${dwi_file}_out_wmfod.nii.gz ${dwi_file}_out_gm ${dwi_file}_out_gmfod.nii.gz ${dwi_file}_out_csf ${dwi_file}_out_csffod.nii.gz -fslgrad ${dwi_file}_mrtrix.bvec ${dwi_file}.bval -force -lmax 6,6,6 -shells 0,1000 -mask ${Patient}_dwi_brainMask.nii.gz
    
    N_tracks=10M
    echo " "
    echo " => TCKGEN "
    tckgen ${dwi_file}_out_wmfod.nii.gz ${Patient}_wm_tracks_${N_tracks}.tck -algorithm iFOD2 -maxlength 600 -seed_gmwmi ${seeding_mask}.nii.gz -cutoff 0.02 -angle 60 -select $N_tracks -backtrack -act ${fivett_file}.nii.gz -force
    tcksift2 ${Patient}_wm_tracks_${N_tracks}.tck ${dwi_file}_out_wmfod.nii.gz ${Patient}_wm_tracks_${N_tracks}_weights.txt -act ${fivett_file}.nii.gz -force


    echo " "
    echo " => Connectome "
    tck2connectome ${Patient}_wm_tracks_${N_tracks}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${N_tracks}.txt -force -zero_diagonal -symmetric -tck_weights_in ${Patient}_wm_tracks_${N_tracks}_weights.txt -assignment_radial_search 10
    tck2connectome ${Patient}_wm_tracks_${N_tracks}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${N_tracks}_density.txt -force -zero_diagonal -symmetric -tck_weights_in ${Patient}_wm_tracks_${N_tracks}_weights.txt -scale_invnodevol -assignment_radial_search 10
    tck2connectome ${Patient}_wm_tracks_${N_tracks}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${N_tracks}_unweighted.txt -force -zero_diagonal -symmetric -assignment_radial_search 10
    tck2connectome ${Patient}_wm_tracks_${N_tracks}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${N_tracks}_density_unweighted.txt -force -zero_diagonal -symmetric -scale_invnodevol -assignment_radial_search 10

done
        

    
