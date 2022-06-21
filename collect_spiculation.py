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


output_dir = 'DATA/LIDC_spiculation'
data_dir = 'DATA/nodule-lidc'
temp_dir = '/tmp'


def task_organize_images(image_file):
    dirname = osp.dirname(image_file)
    dirname = osp.dirname(dirname)
    pid = osp.basename(image_file).split("_")[0]
    print(pid)
    case_output_dir = osp.join(output_dir, pid)
    os.makedirs(case_output_dir, exist_ok=True)
    nodule_mask_files = glob.glob(
        osp.join(dirname, pid, pid+'_CT*-seg-iso-label.nrrd'))
    input_image = sitk.ReadImage(osp.join(image_file))
    for nodule_mask_file in nodule_mask_files:
        print(nodule_mask_file)
        basename = osp.basename(nodule_mask_file)
        input_image_crop_file = nodule_mask_file.replace('-seg-iso-label', '')
        ard_file = nodule_mask_file.replace('label', 'ard')
        ard_surface_file = nodule_mask_file.replace('label', 'ard-surface')
        peaks_file = nodule_mask_file.replace('label', 'peaks-label')
        peaks_surface_file = nodule_mask_file.replace('label', 'peaks-surface')

        mask_image = sitk.ReadImage(nodule_mask_file)
        input_image_crop = sitk.RegionOfInterest(input_image, mask_image.GetSize(
        ), input_image.TransformPhysicalPointToIndex(mask_image.GetOrigin()))
        sitk.WriteImage(input_image_crop, input_image_crop_file, True)

        copylist = [
            {'src': input_image_crop_file, 'dest': osp.join(
                case_output_dir, osp.basename(input_image_crop_file))},
            {'src': nodule_mask_file, 'dest': osp.join(
                case_output_dir, basename)},
            {'src': ard_file, 'dest': osp.join(
                case_output_dir, osp.basename(ard_file))},
            {'src': ard_surface_file, 'dest': osp.join(
                case_output_dir, osp.basename(ard_surface_file))},
            {'src': peaks_file, 'dest': osp.join(
                case_output_dir, osp.basename(peaks_file))},
            {'src': peaks_surface_file, 'dest': osp.join(
                case_output_dir, osp.basename(peaks_surface_file))},
        ]

        for copy in copylist:
            print(copy)
            p = subprocess.Popen(['cp', copy['src'], copy['dest']])
            p.wait()


def task_load_CT_list(ct_data_dir):
    ct_image_list = [osp.join(ct_data_dir, fn)
                     for fn in next(os.walk(ct_data_dir))[2]]
    for ct_image in ct_image_list:
        if ct_image.find('CT') < 0:
            continue

        if ct_image[0] == '.' or ct_image.find('mm') > 0 or \
                ct_image.find('label') > 0 or ct_image.find('csv') > 0 or \
                ct_image.find('txt') > 0 or ct_image.find('levelset') > 0:
            continue

        yield ct_image


if __name__ == "__main__":
    if len(sys.argv) > 1:
        dataset = sys.argv[1]
        print(dataset)
        data_dir = f'DATA/nodule-{dataset.lower()}'
        output_dir = f'DATA/{dataset}_spiculation'
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
        
    for CT_image_file in task_load_CT_list(osp.join(data_dir, 'CT')):
        print(CT_image_file)
        task_organize_images(CT_image_file)
