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
    # Activate conda environment
    module load conda
    conda activate project_env
    Dataset="sample_dataset"
    Parcellation_dir="/workspace/project/resources/parcellation/"
    Template_dir="/workspace/project/resources/template/"
    DATASET_dir="/workspace/project/data/${Dataset}/Nifti/"
fi

Patients=('sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0010')

#'sub-0011' 'sub-0012' 'sub-0013' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0017' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' job_0001 best
#'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0036' 'sub-0037' 'sub-0038' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0043' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047'  job_0002
#'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0053' 'sub-0054' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0066' 'sub-0067' 'sub-0068'  job_0003
#'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0081' 'sub-0082' 'sub-0083' 'sub-0084' 'sub-0085' 'sub-0086' 'sub-0087' 'sub-0088' 'sub-0089' job_0004 best
#'sub-0090' 'sub-0091' 'sub-0092' 'sub-0093' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0100' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0107' 'sub-0108' 'sub-0109' 'sub-0110' 'sub-0111'  job_0005 best
#'sub-0112' 'sub-0113' 'sub-0114' 'sub-0115' 'sub-0116' 'sub-0117' 'sub-0118' 'sub-0119' 'sub-0120' 'sub-0121' 'sub-0122' 'sub-0123' 'sub-0124' 'sub-0125' 'sub-0126' 'sub-0127' 'sub-0128' 'sub-0129' 'sub-0130' 'sub-0131' 'sub-0132' 'sub-0133' 'sub-0134' 'sub-0135' 'sub-0136' job_0006
# 'sub-0137' 'sub-0138' 'sub-0139' 'sub-0140' 'sub-0141' 'sub-0142' 'sub-0143' job_0007
# 'sub-0144' 'sub-0145' 'sub-0146' 'sub-0147' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0153' job_0008 approx 38 min
# 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' job_0009

#'sub-0161' 'sub-0162' job_0010 very slow
#'sub-0163' 'sub-0164' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0172' 'sub-0173' 'sub-0174' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0178' 'sub-0179' 'sub-0180' 'sub-0181' 'sub-0182' 'sub-0183' 'sub-0184' 'sub-0185' 'sub-0186' 'sub-0187' 'sub-0188' 'sub-0189' job_0011 best 
#'sub-0190' 'sub-0191' 'sub-0192' 'sub-0193' 'sub-0194' 'sub-0195' 'sub-0196' 'sub-0197' 'sub-0198' 'sub-0199' 'sub-0200' 'sub-0201' 'sub-0202' 'sub-0203' 'sub-0204' 'sub-0205' 'sub-0206' 'sub-0207' 'sub-0208' 'sub-0209' 'sub-0210' 'sub-0211' 'sub-0212' 'sub-0213' 'sub-0214' 'sub-0215' 'sub-0216' 'sub-0217' 'sub-0218' job_0012 approx 50 min
#'sub-0215' 'sub-0216' 'sub-0217' 'sub-0218' job_0013
#'sub-0219' 'sub-0220' 'sub-0221' 'sub-0222' 'sub-0223' 'sub-0224' 'sub-0225' 'sub-0226' 'sub-0227' 'sub-0228' 'sub-0229' 'sub-0230' 'sub-0231' 'sub-0232' 'sub-0233' 'sub-0234' 'sub-0235' 'sub-0236' 'sub-0237' 'sub-0238' 'sub-0239' 'sub-0240' 'sub-0241' 'sub-0242' 'sub-0243' 'sub-0244' 'sub-0245' 'sub-0246' 'sub-0247' job_0014 approx 50 min
#'sub-0248' 'sub-0249' 'sub-0250' 'sub-0251' 'sub-0252' 'sub-0253' 'sub-0254' 'sub-0255' 'sub-0256' 'sub-0257' 'sub-0258' 'sub-0259' 'sub-0260' 'sub-0261' 'sub-0262' 'sub-0263' 'sub-0264' 'sub-0265' 'sub-0266' 'sub-0267' 'sub-0268' 'sub-0269' 'sub-0270' 'sub-0271' 'sub-0272' 'sub-0273' 'sub-0274' job_0015
#'sub-0275' 'sub-0276' 'sub-0277' 'sub-0278' 'sub-0279' 'sub-0280' 'sub-0281' 'sub-0282' 'sub-0283' 'sub-0284' 'sub-0285' 'sub-0286' 'sub-0287' 'sub-0288' 'sub-0289' 'sub-0290' 'sub-0291' 'sub-0292' 'sub-0293' 'sub-0294' 'sub-0295' 'sub-0296' 'sub-0297' 'sub-0298' 'sub-0299' 'sub-0300' 'sub-0301' 'sub-0302' 'sub-0303' 'sub-0304' 'sub-0305' 'sub-0306' 'sub-0307' 'sub-0308' 'sub-0309' 'sub-0310' job_0016

#'sub-0311' 'sub-0312' 'sub-0313' 'sub-0314' 'sub-0315' 'sub-0316' 'sub-0317' 'sub-0318' 'sub-0319' 'sub-0320' 'sub-0321' 'sub-0322' 'sub-0323' 'sub-0324' 'sub-0325' 'sub-0326' 'sub-0327' 'sub-0328' 'sub-0329' 'sub-0330' 'sub-0331' 'sub-0332' 'sub-0333' 'sub-0334' 'sub-0335' 'sub-0336' 'sub-0337' 'sub-0338' 'sub-0339' 'sub-0340' 'sub-0341' 'sub-0342' 'sub-0343' 'sub-0344' job_0017
#'sub-0345' 'sub-0346' 'sub-0347' 'sub-0348' 'sub-0349' 'sub-0350' 'sub-0351' 'sub-0352' 'sub-0353' 'sub-0354' 'sub-0355' 'sub-0356' 'sub-0357' 'sub-0358' 'sub-0359' 'sub-0360' 'sub-0361' 'sub-0362' 'sub-0363' 'sub-0364' 'sub-0365' 'sub-0366' 'sub-0367' 'sub-0368' 'sub-0369' 'sub-0370' 'sub-0371' 'sub-0372' 'sub-0373' 'sub-0374' 'sub-0375' 'sub-0376' 'sub-0377' 'sub-0378' job_0018
#'sub-0379' 'sub-0380' 'sub-0381' 'sub-0382' 'sub-0383' 'sub-0384' 'sub-0385' 'sub-0386' 'sub-0387' 'sub-0388' 'sub-0389' 'sub-0390' 'sub-0391' 'sub-0392' 'sub-0393' 'sub-0394' 'sub-0395' 'sub-0396' 'sub-0397' 'sub-0398' 'sub-0399' 'sub-0400' 'sub-0401' 'sub-0402' 'sub-0403' 'sub-0404' 'sub-0405' 'sub-0406' 'sub-0407' 'sub-0408' 'sub-0409' 'sub-0410' 'sub-0411' 'sub-0412' job_0019
#'sub-0413' 'sub-0414' 'sub-0415' 'sub-0416' 'sub-0417' 'sub-0418' 'sub-0419' 'sub-0420' 'sub-0421' 'sub-0422' 'sub-0423' 'sub-0424' 'sub-0425' 'sub-0426' 'sub-0427' 'sub-0428' 'sub-0429' 'sub-0430' 'sub-0431' 'sub-0432' 'sub-0433' 'sub-0434' 'sub-0435' 'sub-0436' 'sub-0437' 'sub-0438' 'sub-0439' 'sub-0440' 'sub-0441' 'sub-0442' 'sub-0443' 'sub-0444' 'sub-0445' 'sub-0446' 'sub-0447' 'sub-0448' 'sub-0449' 'sub-0450' 'sub-0451' 'sub-0452' 'sub-0453' 'sub-0454' job_0020

#MCI-
#'sub-0455' 'sub-0456' 'sub-0457' 'sub-0458' 'sub-0459' 'sub-0460' 'sub-0461' 'sub-0462' 'sub-0463' 'sub-0464' 'sub-0465' 'sub-0466' 'sub-0467'  compute-node 4 
#'sub-0468' 'sub-0469' 'sub-0470' 'sub-0471' 'sub-0472' 'sub-0473' 'sub-0474' 'sub-0475' 'sub-0476' 'sub-0477' 'sub-0478' 'sub-0479' 'sub-0480'  compute-node 4 
#'sub-0481' 'sub-0482' 'sub-0147' 'sub-0483' 'sub-0484' 'sub-0485' 'sub-0486' 'sub-0487' 'sub-0488' 'sub-0489' compute-node
#'sub-0490' 'sub-0491' 'sub-0492' 'sub-0493' 'sub-0494' compute-node
#'sub-0495' 'sub-0496' 'sub-0497' 'sub-0289' 'sub-0498' compute-node
#'sub-0499' 'sub-0500' 'sub-0501' 'sub-0502' 'sub-0503' 'sub-0504' 'sub-0505' 'sub-0506' 'sub-0507' 'sub-0508' compute-node
#'sub-0509' 'sub-0510' 'sub-0511' 'sub-0512' 'sub-0513' 'sub-0514'  compute-node
#'sub-0515' 'sub-0516' 'sub-0517' compute-node
#'sub-0518' 'sub-0519'
#'sub-0520' 

#'sub-0521' 'sub-0522' 'sub-0523' 'sub-0524' 'sub-0525' 'sub-0526' 'sub-0527' 'sub-0306' 'sub-0528' 'sub-0529' compute-node job_0021
#'sub-0530' 'sub-0531' 'sub-0532' 'sub-0533' 'sub-0534' 'sub-0535' 'sub-0536' 'sub-0537' 'sub-0538' 'sub-0539' compute-node job_0022
#'sub-0540' 'sub-0541' 'sub-0542' 'sub-0543' 'sub-0544' 'sub-0545' 'sub-0546' 'sub-0547' 'sub-0548' 'sub-0549' compute-node job_0023
#'sub-0550' 'sub-0551' 'sub-0331' 'sub-0552' 'sub-0334' 'sub-0365' 'sub-0553' 'sub-0554' 'sub-0555' 'sub-0556' compute-node job_0024
#'sub-0557' 'sub-0558' 'sub-0559' 'sub-0560' 'sub-0561' 'sub-0562' 'sub-0563' 'sub-0564' 'sub-0565' 'sub-0566' compute-node job_0025
#'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0010' compute-node

for Patient in "${Patients[@]}" 
do
    cd ${DATASET_dir}

    # Initialization 
    #Patient="sub-0043"
    #Patient="sub-0018"
    #Patient_original="${Patients_original[@]}" 
    Cortical_parcellation="Schaefer2018"
    Sub_Cortical_parcellation="Tian_Subcortex_S1"
    Template_Space="MNI152NLin2009cAsym"
    T1_Resolution="1mm"
    N_Cortical_Parcels="400"

    Template_Space2="MIITRA"

    T1_Template_dir=${Template_dir}${Template_Space2}/
    T1_Template=${T1_Template_dir}${Template_Space2}_T1_${T1_Resolution}.nii.gz
    T1_Template_brain=${T1_Template_dir}${Template_Space2}_T1_${T1_Resolution}_brain.nii.gz
    Template_brain_mask=${T1_Template_dir}${Template_Space2}_mask.nii.gz
    T1_subject_brain=${Patient}_T1w_masked
    T1_subject=${Patient}_T1w_reorient.nii.gz
    T1_subject_brain_mask=${Patient}_T1w_brainMask.nii.gz  

    echo "========================================"
    echo "Processing patient:" ${Patient}    

    if [[ "$Cortical_parcellation" == "" ]]  &&  [[ "$Sub_Cortical_parcellation" == "" ]]
    then
        exit "Please provide either a cortical or subcortical parcellation!"
    elif [[ "$Sub_Cortical_parcellation" == "" ]]
    then
        echo "Using only $Cortical_parcellation cortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Cortex/${Cortical_parcellation}_${N_Cortical_Parcels}Parcels_7Networks_order_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks_${Template_Space2}
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks_${Template_Space2}
    elif [[ "Cortical_parcellation" == "" ]]
    then
        echo "Using only $Sub_Cortical_parcellation subcortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Subcortex/${Sub_Cortical_parcellation}_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Sub_Cortical_parcellation}_${Template_Space2}
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Sub_Cortical_parcellation}_${Template_Space2}
    else
        echo "Using $Cortical_parcellation cortical parcellation and $Sub_Cortical_parcellation subcortical parcellation"
        Parcellation=${Parcellation_dir}Schaefer/Cortex-Subcortex/${Cortical_parcellation}_${N_Cortical_Parcels}Parcels_7Networks_order_${Sub_Cortical_parcellation}_${Template_Space}_${T1_Resolution}.nii.gz
        Parcellation_T1_space=${DATASET_dir}${Patient}/anat/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}_${Template_Space2}
        Parcellation_DWI_space=${DATASET_dir}${Patient}/dwi/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}_${Template_Space2}
        Parcellation_fMRI_space=${DATASET_dir}${Patient}/func/${Patient}_space-orig_atlas-${Cortical_parcellation}-${N_Cortical_Parcels}Parcels-7Networks-${Sub_Cortical_parcellation}_${Template_Space2}
    fi

    cd ${DATASET_dir}${Patient}/anat/
    

 
    # Use ANTs for nonlinear registration to template space 
    echo " =>"
    echo "Nonlinear registration to template space ($Template_Space) [ANTS: antsRegistration (Rigid+Affine+NonLinear)]"
    antsRegistration -d 3 --float 0 -o [ ${Patient}_T1to${Template_Space2}_ , ${Patient}_T1to${Template_Space2}.nii.gz] -n Linear -w [ 0.005 , 0.995] -u 0 -r [ $T1_Template_brain, ${T1_subject_brain}.nii.gz, 1] -t Rigid[0.1] -m MI[$T1_Template_brain, ${T1_subject_brain}.nii.gz, 1, 32, Regular, 0.25] -c [ 1000x500x250x100, 1e-7, 10]  -f 8x4x2x1 -s 3x2x1x0vox -t Affine[0.1] -m MI[ $T1_Template_brain, ${T1_subject_brain}.nii.gz, 1, 32, Regular, 0.25] -c [ 1000x500x250x100, 1e-7,10] -f 8x4x2x1 -s 3x2x1x0vox -t SyN[ 0.1, 3, 0] -m CC[ $T1_Template_brain, ${T1_subject_brain}.nii.gz, 1, 4] -c [ 200x200x200x200, 1e-7, 10] -f 8x4x2x1 -s 3x2x1x0vox 
    echo " =>"
    echo "Registering parcellation to subject space [ANTS: antsApplyTransforms] MNI to MIITRA to subject space"
    antsApplyTransforms -d 3 --float 0 -i $Parcellation  -r ${T1_subject_brain}.nii.gz -n NearestNeighbor -t [ ${Patient}_T1to${Template_Space2}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space2}_1InverseWarp.nii.gz -t [ ${T1_Template_dir}Transformation_to_ICBM2009b/Transformation_to_ICBM2009b/MIITRA_to_ICBM2009b_1Affine.txt,1] -t ${T1_Template_dir}Transformation_to_ICBM2009b/Transformation_to_ICBM2009b/MIITRA_to_ICBM2009b_2InverseWarp.nii.gz -t [ ${T1_Template_dir}Transformation_to_ICBM2009b/Transformation_to_ICBM2009b/MIITRA_to_ICBM2009b_3Affine.txt,1] -t ${T1_Template_dir}Transformation_to_ICBM2009b/Transformation_to_ICBM2009b/MIITRA_to_ICBM2009b_4InverseWarp.nii.gz  -o ${Parcellation_T1_space}.nii.gz 
    echo " =>"
    echo "Registering template mask to subject space [ANTS: antsApplyTransforms]"
    antsApplyTransforms -d 3 --float 0 -i $Template_brain_mask  -r ${T1_subject_brain}.nii.gz -n NearestNeighbor  -t [ ${Patient}_T1to${Template_Space2}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space2}_1InverseWarp.nii.gz  -o ${Patient}_T1w_brainMask_tight_${Template_Space2}.nii.gz 
    fslmaths $T1_subject_brain_mask -mul ${Patient}_T1w_brainMask_tight_${Template_Space2}.nii.gz ${Patient}_T1w_brainMask_tight_${Template_Space2}.nii.gz
    echo " =>"
    echo "Registering parcellation to subject space fMRI[FSL: flirt]"
    flirt -in ${Parcellation_T1_space}.nii.gz -ref ${DATASET_dir}${Patient}/func/${Patient}_task-rest_bold_disto_first.nii.gz  -interp nearestneighbour -applyxfm -init ${DATASET_dir}${Patient}/func/${Patient}_T12FMRI_fin.mat -out ${Parcellation_fMRI_space}.nii.gz
    
    ## Register MIITRA priors to subject space and normalize between 0 and 1
    antsApplyTransforms -d 3 --float 0 -i ${T1_Template_dir}${Template_Space2}_csf.nii.gz  -r ${T1_subject_brain}.nii.gz -n Linear -t [ ${Patient}_T1to${Template_Space2}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space2}_1InverseWarp.nii.gz  -o ${T1_subject_brain}_csf_priors_${Template_Space2}.nii.gz
    antsApplyTransforms -d 3 --float 0 -i ${T1_Template_dir}${Template_Space2}_wm.nii.gz  -r ${T1_subject_brain}.nii.gz -n Linear -t [ ${Patient}_T1to${Template_Space2}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space2}_1InverseWarp.nii.gz  -o ${T1_subject_brain}_wm_priors_${Template_Space2}.nii.gz
    antsApplyTransforms -d 3 --float 0 -i ${T1_Template_dir}${Template_Space2}_gm.nii.gz  -r ${T1_subject_brain}.nii.gz -n Linear -t [ ${Patient}_T1to${Template_Space2}_0GenericAffine.mat, 1 ] -t ${Patient}_T1to${Template_Space2}_1InverseWarp.nii.gz  -o ${T1_subject_brain}_gm_priors_${Template_Space2}.nii.gz

    

    cd ${DATASET_dir}${Patient}/dwi/
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse ${DATASET_dir}${Patient}/anat/${T1_subject_brain}_csf_priors_${Template_Space2}.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_csf_priors_${Template_Space2}_coreg.nii.gz -force
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse ${DATASET_dir}${Patient}/anat/${T1_subject_brain}_wm_priors_${Template_Space2}.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_wm_priors_${Template_Space2}_coreg.nii.gz -force
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse ${DATASET_dir}${Patient}/anat/${T1_subject_brain}_gm_priors_${Template_Space2}.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_gm_priors_${Template_Space2}_coreg.nii.gz -force
    mrtransform -linear ${Patient}_DWI2T1_mrtrix.txt -inverse -interp nearest ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg.nii.gz -force
    mrtransform -linear  ${Patient}_DWI2T1_mrtrix.txt -inverse ${Parcellation_T1_space}.nii.gz ${Parcellation_T1_space}_coreg.nii.gz -interp nearest -force

    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_gm_priors_${Template_Space2}_coreg.nii.gz  -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_gm_priors_${Template_Space2}_coreg_axial.nii.gz -R AXIAL
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_wm_priors_${Template_Space2}_coreg.nii.gz  -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_wm_priors_${Template_Space2}_coreg_axial.nii.gz -R AXIAL
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_csf_priors_${Template_Space2}_coreg.nii.gz  -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_csf_priors_${Template_Space2}_coreg_axial.nii.gz -R AXIAL
    animaConvertImage -i ${Parcellation_T1_space}_coreg.nii.gz  -o ${Parcellation_T1_space}_coreg_axial.nii.gz -R AXIAL
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg.nii.gz  -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz -R AXIAL

    
    
    dwi_file=${Patient}_dwi_preprocessed
    T1=${Patient}_T1w_masked_coreg_axial
    fivett_file=${DATASET_dir}${Patient}/anat/${T1}_5tt
    seeding_mask=${DATASET_dir}${Patient}/anat/${T1}_5tt2gmwmi 
    animaConvertImage -i ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial_5tt.nii.gz -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial_5tt_axial.nii.gz -R AXIAL
    fslroi ${fivett_file}_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz 1 1
    fslroi ${fivett_file}_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz 4 1

    cp ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_csf_priors_${Template_Space2}_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_04.nii.gz
    cp ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_wm_priors_${Template_Space2}_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_03.nii.gz
    cp ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_gm_priors_${Template_Space2}_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz

    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz -sub ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz -thr 0 ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz

    CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz 1 1 1
    CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz 1 1 1
    CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_03.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_03.nii.gz 1 1 1
    CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_04.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_04.nii.gz 1 1 1
    CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz 1 1 1
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz -mul ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_02.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz -mul ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_01.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_03.nii.gz -mul ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_03.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_04.nii.gz -mul ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_04.nii.gz
    #fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz -mul ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz
    
    #CopyImageHeaderInformation ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz
    thr=0.35
    thr_suf='35'

    antsAtroposN4.sh -d 3 -a ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_coreg_axial.nii.gz -x ${DATASET_dir}${Patient}/anat/${Patient}_T1w_brainMask_coreg_axial.nii.gz -c 4 -p ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_%02d.nii.gz -w $thr -o ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_

    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation.nii.gz -uthr 1 -bin ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_01.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation.nii.gz -uthr 2 -thr 2 -bin ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_02.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation.nii.gz -uthr 3 -thr 3 -bin ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_03.nii.gz
    fslmaths ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation.nii.gz -uthr 4 -thr 4 -bin ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_04.nii.gz
    cp ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_label_05.nii.gz  ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_05.nii.gz

    dwi_file=${Patient}_dwi_preprocessed
    T1=${Patient}_T1w_masked
    fivett_file=${DATASET_dir}${Patient}/anat/${T1}_${Template_Space2}_coreg_axial_5tt_ants_${thr_suf}
    seeding_mask=${DATASET_dir}${Patient}/anat/${T1}_${Template_Space2}_coreg_axial_5tt2gmwmi_ants_${thr_suf}

    fslmerge -t ${fivett_file}.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_01.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_02.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_03.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_04.nii.gz ${DATASET_dir}${Patient}/anat/${Patient}_T1w_masked_${Template_Space2}_coreg_axial_${thr_suf}_Segmentation_05.nii.gz
    5tt2gmwmi ${fivett_file}.nii.gz ${seeding_mask}.nii.gz -force

    N_tracks=10M
    echo " "
    echo " => TCKGEN "
    tckgen ${dwi_file}_out_wmfod.nii.gz ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck -algorithm iFOD2 -maxlength 600 -seed_gmwmi ${seeding_mask}.nii.gz -cutoff 0.05 -angle 60 -select $N_tracks -backtrack -act ${fivett_file}.nii.gz -force
    tcksift2 ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck ${dwi_file}_out_wmfod.nii.gz ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_weights_${thr_suf}.txt -act ${fivett_file}.nii.gz -force
   
    echo " "
    echo " => Connectome "
    tck2connectome ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${Template_Space2}_${N_tracks}_${thr_suf}.txt -force -zero_diagonal -symmetric -tck_weights_in ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_weights_${thr_suf}.txt
    tck2connectome ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${Template_Space2}_${N_tracks}_density_${thr_suf}.txt -force -zero_diagonal -symmetric -tck_weights_in ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_weights_${thr_suf}.txt -scale_invnodevol
    tck2connectome ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${Template_Space2}_${N_tracks}_unweighted_${thr_suf}.txt -force -zero_diagonal -symmetric
    tck2connectome ${Patient}_${Template_Space2}_wm_tracks_${N_tracks}_${thr_suf}.tck ${Parcellation_T1_space}_coreg_axial.nii.gz connectome_${Template_Space2}_${N_tracks}_density_unweighted_${thr_suf}.txt -force -zero_diagonal -symmetric -scale_invnodevol
done

