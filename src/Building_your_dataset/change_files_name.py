import os
cwd = os.getcwd()

for patient in next(os.walk(cwd))[1]:
    print(patient)
    if 'sub' in patient:
        if os.path.isfile(os.path.join(patient,'fmap',patient+'_fieldmap1.nii.gz')):
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap1.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap1.json'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap3.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap3.json'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap2.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-1_part-phase.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap2.json'),os.path.join(patient,'fmap',patient+'_echo-1_part-phase.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap4.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-2_part-phase.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_fieldmap4.json'),os.path.join(patient,'fmap',patient+'_echo-2_part-phase.json'))
        if os.path.isfile(os.path.join(patient,'fmap',patient+'_echo-1_part-mag_boh.nii.gz')):
            os.rename(os.path.join(patient,'fmap',patient+'_echo-1_part-mag_boh.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-1_part-mag_boh.json'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-2_part-mag_boh.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-2_part-mag_boh.json'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-1_part-phase_boh.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-1_part-phase.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-1_part-phase_boh.json'),os.path.join(patient,'fmap',patient+'_echo-1_part-phase.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-2_part-phase_boh.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-2_part-phase.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_echo-2_part-phase_boh.json'),os.path.join(patient,'fmap',patient+'_echo-2_part-phase.json'))
        if os.path.isfile(os.path.join(patient,'fmap',patient+'_magnitude1.nii.gz')):
            os.rename(os.path.join(patient,'fmap',patient+'_magnitude1.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_magnitude1.json'),os.path.join(patient,'fmap',patient+'_echo-1_part-mag.json'))
            os.rename(os.path.join(patient,'fmap',patient+'_magnitude2.nii.gz'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.nii.gz'))
            os.rename(os.path.join(patient,'fmap',patient+'_magnitude2.json'),os.path.join(patient,'fmap',patient+'_echo-2_part-mag.json'))
