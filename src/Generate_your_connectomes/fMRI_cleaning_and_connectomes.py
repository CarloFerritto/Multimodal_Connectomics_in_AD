import os
import sys
# print(sys.path)  # debug line, can be removed
import nilearn
import numpy as np
from nilearn import plotting
from nilearn.input_data import NiftiMasker
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
import matplotlib.pyplot as plt
import argparse
from io import StringIO

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, required=True, help="Dataset NIFTI file directory")
parser.add_argument('-p', '--patient', type=str, required=True, help='patient id string')
parser.add_argument('-tr', '--TR', type=float, required=True, help='Repetition Time')
args = parser.parse_args()
os.chdir(args.dir)

patient=args.patient
TR=args.TR
func_dir=os.path.join(args.dir,patient,"func")
MIITRA=True
if os.path.isdir(func_dir):
    func_file=os.path.join(func_dir,patient+"_task-rest_bold_norm.nii.gz")
    mask_file=os.path.join(func_dir,patient+"_fMRI_brainMask.nii.gz")
    atlas=os.path.join(func_dir,patient+"_space-orig_atlas-Schaefer2018-400Parcels-7Networks-Tian_Subcortex_S1_MIITRA.nii.gz")
    motion_file = os.path.join(func_dir,patient+"_task-rest_bold_motion_corrected.nii.gz.par")
    if os.path.isfile(func_file) and os.path.isfile(mask_file) and os.path.isfile(atlas) and os.path.isfile(motion_file) and not os.path.isfile(os.path.join(func_dir,"func_connectivity_pearson_correlation_new_MIITRA.txt")):
        nifti_masker = NiftiMasker(mask_img=mask_file)
        time_series=nifti_masker.fit_transform(func_file)
        GS=np.mean(time_series,axis=1)
        GS=np.atleast_2d(GS).T
        confounds = nilearn.image.high_variance_confounds(func_file, n_confounds=5, percentile=2.0, detrend=True, mask_img=None)
        motion = np.loadtxt(motion_file)
        framewise_param= np.diff(motion,axis=0)
        framewise_param[:,:3]=framewise_param[:,:3]*50
        motion_24=np.hstack((motion,np.vstack((np.zeros(6),framewise_param))))
        motion_24=np.hstack((motion_24,motion_24**2)) 
        FD=np.sum(abs(framewise_param),axis=1)
        FD_append=np.atleast_2d(np.append(0,FD)).T
        matrix = np.concatenate((motion_24, confounds,GS,FD_append), axis=1)
        #matrix = np.concatenate((motion, confounds,GS), axis=1)
        FD_new=np.append(0,FD)
        censoring=np.ones((len(FD)+1), dtype=bool)
        for k in range(len(FD_new)):
            if FD_new[k]>0.4:
                if k<1:
                    censoring[:k+3]=False
                elif (k+3)<=len(censoring):
                    censoring[k-1:k+3]=False
                else:
                    censoring[k-1:]=False

        status_switching=np.zeros(len(censoring[censoring==True]))
        h=0
        for k in range(len(censoring)):
            if k>0 and censoring[k]==False and censoring[k-1]==True and h<len(status_switching):
                status_switching[h]=1
            elif censoring[k]==True:
                h+=1

        series_regress=nilearn.signal.clean(time_series, detrend=False,confounds=matrix,low_pass=0.1, high_pass=0.01, t_r=TR,sample_mask=censoring,standardize='zscore_sample')
        rBOLD_regress= nifti_masker.inverse_transform(series_regress)
        rBOLD_regress.to_filename(os.path.join(func_dir,patient+"_task-rest_bold_compcor_new2_MIITRA.nii.gz"))
        masker = NiftiLabelsMasker(labels_img=atlas)
        Roi_time_series=masker.fit_transform(os.path.join(func_dir,patient+"_task-rest_bold_compcor_new2_MIITRA.nii.gz"))

        np.savetxt(os.path.join(func_dir,"time_series_new_MIITRA.txt"),Roi_time_series)
        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([Roi_time_series])[0]
        
        np.fill_diagonal(correlation_matrix, 0)
        
        np.savetxt(os.path.join(func_dir,"func_connectivity_pearson_correlation_new_MIITRA.txt"),correlation_matrix)
        correlation_measure = ConnectivityMeasure(kind='partial correlation')
        correlation_matrix_partial = correlation_measure.fit_transform([Roi_time_series])[0]
        np.fill_diagonal(correlation_matrix_partial, 0)
        np.savetxt(os.path.join(func_dir,"func_connectivity_partial_correlation_new_MIITRA.txt"),correlation_matrix_partial)
        Roi_time_series=Roi_time_series.T
        
        
        f1, (a1, a2, a3, a0) = plt.subplots(4, 1, gridspec_kw={'height_ratios': [1, 1, 1, 4]})
        FD_new=np.append(0,FD)
        a1.plot(FD_new[censoring==True])
        a1.axhline(y=0.4, color='r', linestyle='-')
        #a1.plot(censoring)
        a1.set_title("FD, FD_mean="+str(np.mean(FD))+", FD_max="+str(np.max(FD)))
        a1.autoscale(enable=True, axis='x', tight=True)
        a1.set_xticks([])
        a2.plot(framewise_param[censoring[1:]==True,:])
        a2.set_title("Motion parameters")
        a2.autoscale(enable=True, axis='x', tight=True)
        a2.set_xticks([])
        a3.plot(Roi_time_series[40,:])
        a3.plot(status_switching)
        a3.set_title("Visual network ROI signal after movement regression")
        a3.autoscale(enable=True, axis='x', tight=True)
        a3.set_xticks([])
        a0.imshow(Roi_time_series,aspect='auto',cmap="Greys")
        plt.suptitle("fMRI quality check post censoring: "+patient+", acquisition lenght="+str((np.sum(censoring)-1)*TR)+" seconds")
        f1.set_size_inches(18, 10)
        f1.savefig(os.path.join(func_dir,"fMRI_quality_check_new_MIITRA.png"))

        if os.path.isfile(os.path.join(func_dir,patient+"_task-rest_bold_compcor_new.nii.gz")):
            
            Roi_time_series_old=masker.fit_transform(os.path.join(func_dir,patient+"_task-rest_bold_compcor_new.nii.gz"))
            correlation_measure = ConnectivityMeasure(kind='correlation')
            correlation_matrix_old = correlation_measure.fit_transform([Roi_time_series_old])[0]
            Roi_time_series_old=Roi_time_series_old.T
            np.fill_diagonal(correlation_matrix_old, 0)
            f2, (a4, a5, a6, a7) = plt.subplots(4, 1, gridspec_kw={'height_ratios': [1, 1, 1, 4]})
            a4.plot(FD_new)
            a4.axhline(y=0.4, color='r', linestyle='-')
            a4.plot(censoring)
            a4.set_title("FD, FD_mean="+str(np.mean(FD))+", FD_max="+str(np.max(FD)))
            a4.autoscale(enable=True, axis='x', tight=True)
            a4.set_xticks([])
            a5.plot(framewise_param)
            a5.set_title("Motion parameters")
            a5.autoscale(enable=True, axis='x', tight=True)
            a5.set_xticks([])
            a6.plot(Roi_time_series_old[40,:])
            a6.set_title("Visual network ROI signal after movement regression")
            a6.autoscale(enable=True, axis='x', tight=True)
            a6.set_xticks([])
            a7.imshow(Roi_time_series_old,aspect='auto',cmap="Greys")
            a7.set_title("Carpet plot ROI signals after movement regression")
            a7.set_xlabel("Frame")
            plt.suptitle("fMRI quality check pre censoring: "+patient+", acquisition lenght="+str((len(FD))*TR)+" seconds")
            f2.set_size_inches(18, 10)
            f2.savefig(os.path.join(func_dir,"fMRI_quality_check.png"))
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 20))
            fc_pre_censoring=np.loadtxt(os.path.join(func_dir,"func_connectivity_pearson_correlation.txt"))
            plotting.plot_matrix(fc_pre_censoring, axes=ax[0],title="FC pre censoring",cmap="jet")
            plotting.plot_matrix(correlation_matrix, axes=ax[1],title="FC post censoring",cmap="jet") 
            plotting.plot_matrix(correlation_matrix-fc_pre_censoring,axes=ax[2],title="difference",cmap="jet")
            plt.suptitle("Functional connectivity "+patient)
            fig.set_size_inches(15, 6) 
            fig.savefig(os.path.join(func_dir,"FC_differences.png"))
        if MIITRA:
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 20))
            fc_MNI=np.loadtxt(os.path.join(func_dir,"func_connectivity_pearson_correlation_new.txt"))
            plotting.plot_matrix(fc_MNI, axes=ax[0],title="FC MNI",cmap="jet")
            plotting.plot_matrix(correlation_matrix, axes=ax[1],title="FC MIITRA",cmap="jet") 
            plotting.plot_matrix(correlation_matrix-fc_MNI,axes=ax[2],title="difference",cmap="jet")
            plt.suptitle("Functional connectivity "+patient)
            fig.set_size_inches(15, 6) 
            fig.savefig(os.path.join(func_dir,"FC_differences_MIITRA.png"))


