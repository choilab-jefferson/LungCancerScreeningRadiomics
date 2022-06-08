import os
import subprocess

import os.path as osp
import numpy as np
import SimpleITK as sitk

from utils import nodule_info


def image_crop_by_mask(image_file, output_file, mask_file):
    input_image = sitk.ReadImage(image_file)
    mask_image = sitk.ReadImage(mask_file)
    input_image_crop = sitk.RegionOfInterest(input_image, mask_image.GetSize(), input_image.TransformPhysicalPointToIndex(mask_image.GetOrigin()))
    sitk.WriteImage(input_image_crop, output_file, True)


def image_resample(input_file, output_file, iso_size="1"):
    image_resample = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools python PythonTools/image_resample.py".split(" ")
    if type(input_file) == list:
        input_file = input_file[0]
    
    p = subprocess.Popen(image_resample+ [osp.abspath(input_file), osp.abspath(output_file), iso_size, iso_size, iso_size])
    #p = subprocess.Popen(["touch", output_file])
    p.wait()


def extract_nodule_labels(input_file, output_file):
    print(input_file)
    print(output_file)
    dirname = osp.dirname(input_file)
    names = osp.basename(input_file).split("_")
    basename = names[0]
    nidx = int(names[2].replace('.nodule', '')) - 1

    input_nodule_info = osp.join(dirname, basename + '_all.csv')
    nodules = nodule_info.Nodules(input_nodule_info).nodules

    nodule = nodules.loc[nidx, ]
    bb_idx = range(nodule.index.get_loc('BoundingBox_idx_1'),
                   nodule.index.get_loc('BoundingBox_idx_6') + 1)
    bbox = np.asarray(nodule.iloc[bb_idx]).astype(np.uint32)
    bbox_size = bbox[3:6].copy()
    bbox[3:6] = bbox[0:3] + bbox[3:6]
    print(basename, nidx + 1, bbox)

    new_mask_file = ['', '', '', '']
    for mi in range(4):
        sid = 'Phy' + str(mi + 1)
        mask_file = osp.join(dirname, basename + '_CT_' + sid + '-label.nrrd')
        new_mask_file[mi] = osp.join(
            dirname, basename + '_CT_' + str(nidx + 1) + '-' + sid + '-label.nrrd')
        if osp.isfile(mask_file) == False:
            with open(new_mask_file[mi], "w"):
                pass
            print(basename, nidx + 1, sid, bbox, 'no1')
        else:
            mask_image = sitk.ReadImage(mask_file)
            mask_size = mask_image.GetSize()
            mask_spacing = np.asarray(mask_image.GetSpacing())
            bbox1 = bbox.copy()
            bbox1[5] = mask_size[2] - bbox[2]
            bbox1[2] = mask_size[2] - bbox[5]
            pad = (4 / mask_spacing).astype(int)
            crop_o = np.asarray(
                [np.asarray(bbox1[0:3]), mask_size - bbox1[3:6]])
            crop = np.asarray([np.asarray(bbox1[0:3]) - pad,
                               mask_size - np.round(bbox1[3:6] + pad)])
            crop[np.asarray(crop) < 0] = 0
            print(bbox, bbox_size, mask_size, crop_o, crop)
            center_image = sitk.Crop(
                mask_image, crop_o[0].tolist(), crop_o[1].tolist())
            mask_image = sitk.Crop(
                mask_image, crop[0].tolist(), crop[1].tolist())

            filter = sitk.LabelShapeStatisticsImageFilter()
            filter.Execute(center_image)
            if filter.GetNumberOfLabels() == 0:
                with open(new_mask_file[mi], "w"):
                    pass
                print(basename, nidx + 1, sid, bbox1, 'no2')
            else:
                print(basename, nidx + 1, sid, bbox1)
                sitk.WriteImage(mask_image, new_mask_file[mi], True)


def staple_comparison(input_file, output_file):
    staple_comparison = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools STAPLEComparison".split(" ")
    STAPLE_to_labelmap = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools python PythonTools/STAPLE_to_labelmap.py".split(" ")
    input_file = [osp.abspath(i) for i in input_file]
    output_file = [osp.abspath(o) for o in output_file]
    print(input_file)
    print(output_file)
    dirname = osp.dirname(input_file[0])
    basename = osp.basename(input_file[0]).split("_")[0]

    for mi in range(4):
        if osp.getsize(input_file[mi]) == 0:
            input_file[mi] = ''

    new_mask_file = np.delete(input_file, np.where(
        np.asarray(input_file) == '')).tolist()

    p = subprocess.Popen(staple_comparison + [output_file[1], output_file[2], '0'] + new_mask_file)
    p.wait()
    p = subprocess.Popen(STAPLE_to_labelmap + [output_file[1], output_file[0], '0.5'])
    p.wait()

