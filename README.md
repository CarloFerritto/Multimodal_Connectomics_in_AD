# Multimodal Connectomics in AD(Ferritto et al.)

Brain alterations in Alzheimer's disease are often studied through either structural or functional MRI, even though these two modalities provide different and complementary information. This study shows that combining them helps capture disease-related network changes more effectively than considering either one alone. Importantly, the advantage of multimodal models is associated with brain regions and networks that are consistent with current knowledge of Alzheimer's disease. It also suggests that multimodal connectivity models are biologically interpretable, helping relate predictive results to specific brain regions and networks rather than to purely statistical effects.
  
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

See the pre-print on [HAL website](https://hal.science/hal-05649878).

# Contents overview
In this repository you can find instructions and scripts to reproduce our results. <br>
There's no mean to reproduce the exact same results, as we can't provide the subject ID's and visits used in this study. <br>
Nevertheless, we provide specific instruction on how to select the set of subjects.



## Reproducing full analysis
### Building your Dataset 
In this section we provide istruction on how to build your dataset. All the Excel file can be downloaded [here](https://ida.loni.usc.edu).<br>
Before dowloading the images we suggest to:
   * filter the patients through the advance search tool and select: ADNI3 as cohort.
   * filter the patients through the "Mayo(Jack Lab) - ADNI 3 MRI QC" file, selecting only the patients the have the T1w, DWI, fMRI and field map (optional) which passed the quality check.
   * filter the patients through "Diagnosis" and "AMYLOID STATUS" (in the UCBERKELEY_AMY file).

Once you have all the possible subjects we suggest to build an excel file similar to [dataset.xlsx](data/dataset.xlsx). <br>
Based on the resulting Excel file, only one visit per patient should be selected, and only visits with all values available. <br>


### Images Download and organization
Images can be download  [here](https://ida.loni.usc.edu).<br>
Once you download the images, one should:
 * convert the images in Nifti format ( we suggest to use the [heudiconv tool](https://github.com/nipy/heudiconv) with the [heuristic.py](src/Building_your_dataset/heuristic.py) and [heuristic2.py](src/Building_your_dataset/heuristic2.py) scripts and use the [change_files_name.py](src/Building_your_dataset/change_files_name.py) script to change images string)
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
      ├─ func/
      │  ├─ sub-003_task-rest_bold.nii.gz
      │  └─ sub-003_task-rest_bold.json
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
   │  └─ Transformation_to_ICBM2009b/
   │     ├─ MIITRA_to_ICBM2009b_1Affine.txt
   │     ├─ MIITRA_to_ICBM2009b_2Warp.txt
   │     ├─ MIITRA_to_ICBM2009b_2InverseWarp.txt
   │     ├─ MIITRA_to_ICBM2009b_3Affine.txt
   │     ├─ MIITRA_to_ICBM2009b_4Warp.txt
   │     └─ MIITRA_to_ICBM2009b_4InverseWarp.txt
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

To run the preprocessing use the [pre_processing_pipeline.sh](src/Preprocess_your_images/pre_processing_pipeline.sh) script (be sure to have the [metadata_handler.py](src/Preprocess_your_images/metadata_handler.py) and [flip_bvec.py](src/Preprocess_your_images/flip_bvec.py) script in tha same folder).

### Generate your connectomes 
To generate the structural connectomes use the [tractography_pipeline.sh](src/Generate_your_connectomes/tractography_pipeline.sh) script, while for the functional connectomes use the [fMRI_cleaning_and_connectomes.py](src/Generate_your_connectomes/fMRI_cleaning_and_connectomes.py).

### Harmonize your connectomes
To harmonize your connectomes, first select the calibration set using the [select_calibration_set.py](src/Harmonize_your_connectomes/select_calibration_set.py) script and then harmonize using the [harmonize_connectomes.py](src/Harmonize_your_connectomes/select_calibration_set.py) script.
### Analyze your connectome
To compute the nodal connectome metrics use the [Metrics_extraction.py](src/Analyze_your_connectomes/Metrics_extraction.py) script and the [environment_preprocessing_and_metrics.yaml](src/environment_preprocessing_and_metrics.yml) environment. <br>
To run classification and regression analysis use the [metrics_extraction.py](src/Analyze_your_connectomes/metrics_extraction.py) and [regression_analysis.py](src/Analyze_your_connectomes/regression_analysis.py) scripts and the [environment_analysis.yaml](src/environment_analysis.yaml) environment.
To reproduce the results and the plots use the [results_classification.ipynb](src/Analyze_your_connectomes/results_classification.ipynb), [results_regression.py](src/Analyze_your_connectomes/results_regression.ipynb)  and [results_overlap_analysis_classification_regression.py](src/Analyze_your_connectomes/results_overlap_analysis_classification_regression.py) scripts.





