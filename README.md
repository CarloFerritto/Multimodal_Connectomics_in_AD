# Multimodal Connectomics in Alzheimer's Disease  (Ferritto et al.)

In this study, we investigated how preprocessing choices shape late life structural connectomics in cognitively healty subjects.<br>
We held tractography and parcellation constant while varying two upstream factors: the reference template and the tissue segmentation strategy that drives anatomically constrained tractography.
  
## Table of contents
   * [How to cite?](#how-to-cite)
   * [Contents overview](#contents-overview)
   * [Reproducing full analysis](#reproducing-full-analysis)
      * [Building your dataset](#building-your-dataset)
      * [Images download and organization](#images-download-and-organization)
      * [Preprocess your images](#preprocess-your-images)
      * [Generate your connectomes](#generete-your-connectomes)
      * [Analyze your connectomes](#analyze-your-connectome)

## How to cite?


# Contents overview
In this repository you can find instructions and scripts to reproduce our results. <br>
There's no mean to reproduce the exact same results, as we can't provide the subject ID's and visits used in this study. <br>
Nevertheless, we provide specific instruction on how to select the set of subjects.



## Reproducing full analysis
### Building your Dataset 
I

### Images Download and organization
Images can be download  [here](https://ida.loni.usc.edu).<br>
Once you download the images, one should:
 * convert the images in Nifti format ( we suggest to use the [heudiconv tool](https://github.com/nipy/heudiconv) with the [heuristic.py](src/Building_your_dataset/heuristic.py) script and use the [change_files_name.py](src/Building_your_dataset/change_files_name.py) script to change images string)
 * organize the folder as follow 
```text
.
└─ Nifti/
   ├─ sub-001/
   ├─ sub-002/
   └─ sub-003/
      ├─ anat/
      │  ├─ sub-003_T1w.nii.gz
      │  └─ sub-003_T1w.json
      ├─ dwi/
      │  ├─ sub-003_dwi.nii.gz
      │  ├─ sub-003_dwi.bval
      │  ├─ sub-003_dwi.bvec
      │  └─ sub-003_dwi.json
      └─ fmap/
         ├─ sub-003_phasediff.nii.gz
         ├─ sub-003_phasediff.json
         ├─ sub-003_echo-1_part-mag.nii.gz
         ├─ sub-003_echo-1_part-mag.json
         ├─ sub-003_echo-2_part-mag.nii.gz
         └─ sub-003_echo-2_part-mag.json
```

### Preprocess your images 
To preprocess the images you need to dowload [Anima](https://anima.readthedocs.io/en/latest/), [ANTs](https://github.com/ANTsX/ANTs) 2.6.0.dev1-gb775a15, [FSL](https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FslInstallation.html) 6.0.7.17, and use the [environment_preprocessing_and_metrics.yml](src/environment_preprocessing_and_metrics.yml) environment. Be sure to download the [MIITRA](https://www.nitrc.org/frs/?group_id=1407) and [MNI152NLin2009cAsym](https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009) templates, and the [400 Schaefer+ S1 Tian](https://github.com/yetianmed/subcortex) parcellation.<br>
One should have:
 * one folder containing the templates
```text
.
└─ TEMPLATE/
   ├─ MIITRA/                       
   │  ├─ MIITRA_T1_1mm.nii.gz
   │  ├─ MIITRA_T1_1mm_brain.nii.gz
   │  ├─ MIITRA_mask.nii.gz
   │  ├─ MIITRA_gm.nii.gz
   |  ├─ MIITRA_csf.nii.gz
   │  └─ MIITRA_wm.nii.gz
   └─ MNI152NLin2009cAsym/          
      ├─ MNI152NLin2009cAsym_T1_1mm.nii.gz
      ├─ MNI152NLin2009cAsym_T1_1mm_brain.nii.gz
      ├─ MNI152NLin2009cAsym_mask.nii.gz
      ├─ MNI152NLin2009cAsym_gm.nii.gz
      ├─ MNI152NLin2009cAsym_csf.nii.gz
      └─ MNI152NLin2009cAsym_wm.nii.gz
```
 * one folder containing the parcellation
```text
.
└─ PARCELLATION/
   ├─ Schaefer/
   └─ Cortex-Subcortex/             
      ├─ Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S1_MNI152NLin2009cAsym_1mm.nii.gz
      └─ Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S1_MNI152_label.txt
```



### Generate your connectomes 


### Analyze your connectome






