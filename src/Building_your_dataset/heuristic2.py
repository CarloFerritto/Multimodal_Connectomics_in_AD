
import os
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    T1w = create_key('sub-{subject}/anat/sub-{subject}_T1w')
    fMRI = create_key('sub-{subject}/func/sub-{subject}_task-rest_bold')
    dwi = create_key('sub-{subject}/dwi/sub-{subject}_dwi')
    fmap_mag =  create_key('sub-{subject}/fmap/sub-{subject}_magnitude')
    fmap_phase = create_key('sub-{subject}/fmap/sub-{subject}_phasediff')
    fmap_numb = create_key('sub-{subject}/fmap/sub-{subject}_boh')
    flair=create_key('sub-{subject}/anat/sub-{subject}_FLAIR')
    info = {T1w: [], fMRI: [], dwi: [], fmap_mag: [], fmap_phase: [], fmap_numb: [],flair: [],}
    
    for idx, s in enumerate(seqinfo):
        if ('Accel' in s.series_description) or ('accel' in s.series_description):
            info[T1w].append(s.series_id)
        if ('Axial' in s.series_description) and ('MRI' in s.series_description):
            info[fMRI].append(s.series_id)            
        if ('DTI' in s.series_description):
            info[dwi].append(s.series_id)
        if ((s.dim3 == 108) or (s.dim3 == 90)) and ('Field' in s.series_description) and not('Axial' in s.series_description):
            info[fmap_mag].append(s.series_id)
        if ((s.dim3 == 54) or (s.dim3 == 45)) and ('Field' in s.series_description) and not('Axial' in s.series_description):
            info[fmap_phase].append(s.series_id)
        if ('Axial Field Mapping' in s.series_description):
            info[fmap_numb].append(s.series_id)
        if ('FLAIR' in s.series_description):
            info[flair].append(s.series_id)

    return info

