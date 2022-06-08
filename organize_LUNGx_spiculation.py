#!/usr/bin/env python3
import os
import sys
import glob
import time
import subprocess

import os.path as osp
import scipy as sp
import numpy as np
import SimpleITK as sitk


output_dir = 'DATA/lungx/LUNGx_spiculation'
data_dir = 'DATA/lungx/LUNGx_CT'
temp_dir = 'tmp'


def task_organize_images(image_file):
    dirname = osp.dirname(image_file)
    dirname = osp.dirname(dirname)
    pid = osp.basename(image_file).split("_")[0]
    print(pid)
    case_output_dir = osp.join(output_dir, pid)
    os.makedirs(case_output_dir, exist_ok=True)
    nodule_mask_files = glob.glob(osp.join(dirname,'radiomics',pid+'_CT*-seg-label.nrrd'))
    #input_image = sitk.ReadImage(osp.join(dirname,'LUNGx_CT',pid+'_CT.nrrd'))
    for nodule_mask_file in nodule_mask_files:
        print(nodule_mask_file)
        basename = osp.basename(nodule_mask_file)
        input_image_crop_file = nodule_mask_file.replace('-seg-label', '')
        
        #mask_image = sitk.ReadImage(nodule_mask_file)
        #input_image_crop = sitk.RegionOfInterest(input_image, mask_image.GetSize(), input_image.TransformPhysicalPointToIndex(mask_image.GetOrigin()))
        #sitk.WriteImage(input_image_crop, input_image_crop_file, True)

        copylist = [
            {'src': input_image_crop_file, 'dest': osp.join(case_output_dir, basename.replace("seg-label","gc+cip"))},
            {'src': nodule_mask_file, 'dest': osp.join(case_output_dir, basename.replace("seg","gc+cip"))},
        ]

        for copy in copylist:
            print(copy)
            p = subprocess.Popen(['cp', copy['src'], copy['dest']])
            p.wait()
    
    
def task_load_dicom_list(data_dir):
    org_patients_list = [osp.join(data_dir, fn) for fn in next(os.walk(data_dir))[2]]
    for patient in org_patients_list:
        pid = osp.basename(patient)
        #print(pid)
        if pid.find('CT') < 0:
            continue
        
        if pid[0] == '.' or pid.find('mm') > 0 or pid.find('label') > 0 or pid.find('csv') > 0 or pid.find('txt') > 0 or pid.find('levelset') > 0:
            continue

        pt_nrrd_path = osp.join(data_dir, pid)
        yield pt_nrrd_path


for CT_image_file in task_load_dicom_list(data_dir):
    print(CT_image_file)
    task_organize_images(CT_image_file)