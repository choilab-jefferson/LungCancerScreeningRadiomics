#!/usr/bin/env python

import os
import sys
import glob
import time
import subprocess

import os.path as osp
import scipy as sp
import numpy as np
import pandas as pd
import SimpleITK as sitk

from multiprocessing import freeze_support
from ruffus import *
import ruffus.cmdline as cmdline
from utils import metadata

import tasks
from tasks.segmentation import _image_resample, _cip_segmentation, _compare_masks

##########################################################################
# Set your environmental parameters
##########################################################################
output_dir = 'output'
data_dir = 'DATA'
dicom_path = 'DATA/LUNGx'
iso_size = '1'

experiment_set = 'nodule-lungx'
meta = metadata('metadata/LUNGx.csv')
##########################################################################


def retrive_dicom_list():
    # retrieve dicom dirs and preparation of the input image sets
    patient_path_list = [os.path.join(dicom_path, fn) for fn in next(os.walk(dicom_path))[1]]

    for pt_dicom_path in patient_path_list:
        pid = os.path.basename(pt_dicom_path)

        pt = meta.getPatient(pid)
        if len(pt) == 0:
            continue

        yield pt_dicom_path


def task_segment_nodule_growcut(input_file, output_files, pid, output_prefix):
    nodule_segmentation = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools NoduleSegmentation".split(" ")
    input_file = osp.abspath(input_file)
    output_prefix = osp.abspath(output_prefix)
    output_files = [osp.abspath(o) for o in output_files]
    pid = pid.strip()
    pt = meta.getPatient(pid)
    print(pid, pt, output_files, output_prefix)
    os.makedirs(osp.dirname(output_prefix), exist_ok=True)

    for idx in pt:
        pt_row = pt[idx]
        x = pt_row['X']
        y = pt_row['Y']
        slice_number = pt_row['Z']
        large_diameter = pt_row['LD']
        short_diameter = pt_row['PD']
        print(idx, slice_number, large_diameter, short_diameter)

        iso_output_file = output_prefix + str(idx) + "-iso-label.nrrd"
        p = subprocess.Popen(nodule_segmentation + [input_file, str(x), str(y), str(slice_number), str(large_diameter), str(short_diameter), str(iso_output_file)])
        p.wait()

        mask_image = sitk.ReadImage(iso_output_file)
        input_image = sitk.ReadImage(iso_output_file+'in.nrrd')
        mask_image = sitk.VotingBinaryIterativeHoleFilling(mask_image)
        mask_image = _image_resample(mask_image, input_image, True) # to make it in the same space
        output_file = output_prefix + str(idx) + "-label.nrrd"
        sitk.WriteImage(mask_image, output_file, True)


def task_segment_nodule(in_mask_file, output_files, image_file):
    nodule_segmentation_cip = f"docker run --volume {os.getcwd()}:{os.getcwd()} acilbwh/chestimagingplatform GenerateLesionSegmentation".split(" ")
    mask_file, levelset_file = output_files[0], output_files[1]

    input_image = sitk.ReadImage(in_mask_file.replace('-label','-iso-label.nrrdin'))
    mask_image = sitk.ReadImage(in_mask_file)

    filter = sitk.LabelShapeStatisticsImageFilter()
    filter.ComputeFeretDiameterOn()
    filter.Execute(mask_image)
    center = filter.GetCentroid(1)
    radius = filter.GetFeretDiameter(1) / 2
    print(center)
    print(radius)

    dirname = osp.dirname(in_mask_file)
    names = osp.basename(in_mask_file).split("_")
    pid = names[0]

    np_mask_image = sitk.GetArrayFromImage(mask_image) > 0
    np_input_image = sitk.GetArrayFromImage(input_image)

    print(mask_image.GetSize())
    print(input_image.GetSize())

    ## chest image platform (CIP) segmenation
    partSolid = '--echo'
    # if nodule.loc['Texture'] != 5:
    #    partSolid = '--partSolid'
    seeds = list()
    seeds.append('--seeds')
    seeds.append('%f,%f,%f' % center)
    try:
        #points = np.matrix(np.where(np.bitwise_and(np_staple_image > 0.9, np_input_image > np.percentile(np_input_image[np_mask_image], 95))))
        points = np.matrix(np.where(np.bitwise_and(np_mask_image, np_input_image > np.percentile(np_input_image[np_mask_image], 90))))
        seed_points = points[:, np.random.choice(range(0, points.shape[1]), size=min(points.shape[1], 100), replace=False)]

        for seed in seed_points.transpose().tolist():
            seed_phy = mask_image.TransformIndexToPhysicalPoint(
                (seed[2], seed[1], seed[0]))
            seeds.append('--seeds')
            seeds.append('%f,%f,%f' % seed_phy)
            # print(seed_phy)
    except:
        print(seeds)

    print(levelset_file)
    for _ in range(5):
        try:
            _cip_segmentation(image_file, levelset_file, radius, seeds, partSolid)
            levelset_image = sitk.ReadImage(levelset_file)
            break
        except:
            time.sleep(1)
            continue

    sys.stdout.flush()


    ## grow cut segmenation
    gc_mask_image = sitk.ReadImage(in_mask_file.replace('-label','-iso-label'))
    """
    gc_mask_file = mask_file.replace('seg', 'growcut')
    center_idx = input_image.TransformPhysicalPointToIndex(center)
    print(center_idx)
    p = subprocess.Popen([nodule_segmentation, image_file, str(center_idx[0]),
                          str(center_idx[1]), str(-center_idx[2]), str(radius*2+5), str(radius), gc_mask_file])
    p.wait()

    gc_mask_image = sitk.ReadImage(gc_mask_file)
    """

    gc_mask_image = sitk.Cast(sitk.BinaryThreshold(gc_mask_image, 1, 1000), sitk.sitkInt32)
    gc_mask_image = _image_resample(gc_mask_image, levelset_image, True) # to make it in the same space
    gc_mask_image = sitk.ConstantPad(gc_mask_image, (100, 100, 100), (100, 100, 100))
    gc_mask_image = sitk.RegionOfInterest(gc_mask_image, levelset_image.GetSize(),
        gc_mask_image.TransformPhysicalPointToIndex(levelset_image.GetOrigin()))
    gc_mask_image.SetOrigin(levelset_image.GetOrigin())



    ## make physician's contour to isotropic voxel
    mask_image = sitk.Cast(sitk.BinaryThreshold(mask_image, 1, 1000), sitk.sitkInt32)
    iso_mask_image = _image_resample(mask_image, levelset_image, True)
    iso_mask_image = sitk.ConstantPad(iso_mask_image, (100, 100, 100), (100, 100, 100))
    iso_mask_image = sitk.RegionOfInterest(iso_mask_image, levelset_image.GetSize(),
        iso_mask_image.TransformPhysicalPointToIndex(levelset_image.GetOrigin()))
    iso_mask_image.SetOrigin(iso_mask_image.GetOrigin())
    sitk.WriteImage(iso_mask_image, mask_file.replace('seg', 'all-iso'), True)


    ## post processing
    cip_mask_image = sitk.BinaryThreshold(levelset_image, 0)
    cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    # cip_mask_image = sitk.BinaryClosingByReconstruction(
    #     cip_mask_image, 1, sitk.sitkBall, 1, True)
    # cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    cip_mask_image = sitk.Cast(cip_mask_image, sitk.sitkInt32)

    gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)
    # gc_mask_image = sitk.BinaryClosingByReconstruction(
    #     gc_mask_image, 1, sitk.sitkBall, 1, True)
    # gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)


    fsegmesure = open(mask_file.replace('-label.nrrd', '.txt'),'w')


    union_image = sitk.Or(cip_mask_image, gc_mask_image)
    intersection_image = sitk.And(cip_mask_image, gc_mask_image)

    ## select bigger and smaller segmentation from growcut and cip
    U = (sitk.GetArrayFromImage(union_image)>0).sum()
    I = (sitk.GetArrayFromImage(intersection_image)>0).sum()
    G = (sitk.GetArrayFromImage(gc_mask_image)>0).sum()
    C = (sitk.GetArrayFromImage(cip_mask_image)>0).sum()
    G_ = G-I
    C_ = C-I

    print(U,I,G,C,I/U,G/U,C/U,G_,C_)
    print(U,I,G,C,I/U,G/U,C/U,G_,C_, file=fsegmesure)

    if G_ > C_:
        bigger_mask_image,smaller_mask_image = gc_mask_image, cip_mask_image
    else:
        bigger_mask_image,smaller_mask_image = cip_mask_image, gc_mask_image
    if (G_<=0 or C_<=0) and I/U<0.5:
        smaller_mask_image = bigger_mask_image

    sitk.WriteImage(gc_mask_image, mask_file.replace('seg', 'gc-iso'), True)
    sitk.WriteImage(cip_mask_image, mask_file.replace('seg', 'cip-iso'), True)


    ## compare with physician's contour
    aniso_gc_mask_image = _image_resample(gc_mask_image, mask_image, True)
    aniso_cip_mask_image = _image_resample(cip_mask_image, mask_image, True)

    mask_image = sitk.ConstantPad(mask_image, (100, 100, 100), (100, 100, 100))
    mask_image = sitk.RegionOfInterest(mask_image, aniso_gc_mask_image.GetSize(),
        mask_image.TransformPhysicalPointToIndex(aniso_gc_mask_image.GetOrigin()))
    mask_image.SetOrigin(aniso_gc_mask_image.GetOrigin())

    sitk.WriteImage(aniso_gc_mask_image, mask_file.replace('seg', 'gc'), True)
    sitk.WriteImage(aniso_cip_mask_image, mask_file.replace('seg', 'cip'), True)

    gc_compare = _compare_masks(mask_image, aniso_gc_mask_image)
    cip_compare = _compare_masks(mask_image, aniso_cip_mask_image)
    print("Growcut", file=fsegmesure)
    print(gc_compare, file=fsegmesure)
    print("CIP", file=fsegmesure)
    print(cip_compare, file=fsegmesure)


    ## find attached area
    subt_image = sitk.Subtract(union_image, intersection_image)
    subt_image = sitk.BinaryOpeningByReconstruction(subt_image, [1,1,1], sitk.sitkBall, 1, True)

    mask_labelmap = sitk.BinaryImageToLabelMap(subt_image)
    mask_label = sitk.LabelMapToLabel(mask_labelmap)
    filter.Execute(mask_label)
    max_size = 0
    s = 0
    for l in filter.GetLabels():
        size = filter.GetNumberOfPixels(l)
        print("Label: {0}, Size: {1}".format(l, size))
        if size > max_size:
            s = l
            max_size = size

    attched_image = sitk.LabelMapMask(mask_labelmap, union_image, s)
    sitk.WriteImage(attched_image, mask_file.replace(
        'seg', 'attached'), True)


    ## remove attached area from the segmenation
    temp_image = sitk.LabelMapMask(mask_labelmap, smaller_mask_image, s, negated=True)
    mask_labelmap = sitk.BinaryImageToLabelMap(temp_image)
    mask_label = sitk.LabelMapToLabel(mask_labelmap)
    filter.Execute(mask_label)
    max_size = 0
    s = 0
    for l in filter.GetLabels():
        size = filter.GetNumberOfPixels(l)
        print("Label: {0}, Size: {1}".format(l, size))
        if size > max_size:
            s = l
            max_size = size
    out_mask_image = sitk.LabelMapMask(mask_labelmap, temp_image, s)
    sitk.WriteImage(out_mask_image, mask_file.replace('seg', 'seg-iso'), True)

    aniso_out_mask_image = _image_resample(out_mask_image, mask_image, True)
    aniso_out_mask_image.SetOrigin(mask_image.GetOrigin())
    sitk.WriteImage(aniso_out_mask_image, mask_file, True)

    seg_compare = _compare_masks(mask_image, aniso_out_mask_image)
    print("Segmentation", file=fsegmesure)
    print(seg_compare, file=fsegmesure)


    fsegmesure.close()
    # attched_image = sitk.And(sitk.BinaryDilate(
    #    output_image), sitk.BinaryDilate(attched_obj_image))
    #sitk.WriteImage(attched_image, mask_file.replace('seg', 'attached'), True)


def make_pipeline_lungx(experiment_set):
    pipeline_name = experiment_set
    data_path = osp.join(data_dir, experiment_set)
    output_path = osp.join(output_dir, experiment_set)
    feature_list_path = osp.join(output_dir, "feature-list_" + experiment_set + ".csv")
    pipeline = Pipeline(pipeline_name)
    
    dicom_list = list(retrive_dicom_list())

    pipeline.transform(name="task_dicom_to_nrrd_convert",
                       task_func=tasks.dicom_to_nrrd_convert,
                       input=dicom_list,
                       filter=formatter(),
                       output=data_path + "/CT/{basename[0]}_CT.nrrd")\
        .follows(mkdir(output_path)) \
        .follows(mkdir(data_path)) \
        .follows(mkdir(data_path+"/CT"))


    pipeline.subdivide(name="task_segment_nodule_growcut",
                       task_func=task_segment_nodule_growcut,
                       input=output_from("task_dicom_to_nrrd_convert"),
                       filter=formatter(experiment_set + r"/CT/(?P<pid>[^/]+)_CT.nrrd"),
                       output=data_path + "/{pid[0]}/{basename[0]}_*[0-9]-label.nrrd",
                       # after '_' will be index of nodule which is some digit
                       extras=["{pid[0]}", data_path + "/{pid[0]}/{basename[0]}_"]) \
                           

    pipeline.transform(name="task_segment_nodule",
                       task_func=task_segment_nodule,
                       input=output_from("task_segment_nodule_growcut"),
                       filter=formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-label.nrrd"),
                       output=[data_path + "/{pid[0]}/{image_name[0]}_{nid[0]}-seg-label.nrrd",
                              data_path + "/{pid[0]}/{image_name[0]}_{nid[0]}-levelset-label.nrrd"],
                       extras=[data_path + "/CT/{pid[0]}_CT.nrrd"])

    pipeline.merge(name="task_segmentation_organization",
                   task_func=tasks.segmentation_organization,
                   input=output_from("task_segment_nodule"),
                   output=output_dir + "/Segmentation_" + experiment_set + ".csv")

    pipeline.transform(name="task_image_resample",
                       task_func=tasks.image_resample,
                       input=output_from("task_segment_nodule_growcut"),
                       filter=formatter(),
                       output="{path[0]}/{basename[0]}-" + iso_size + "mm.nrrd")

    pipeline.transform(name="task_feature_extraction",
                       task_func=tasks.feature_extraction,
                       input=output_from("task_image_resample"),
                       filter=formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-label-" + iso_size + "mm.nrrd"),
                       output=output_path + "/{image_name[0]}_{nid[0]}-" + iso_size + "mm.txt",
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"]) \
        .follows(mkdir(output_dir), mkdir(output_path))

    pipeline.transform(name="task_image_resample_seg",
                       task_func=tasks.image_resample,
                       input=output_from("task_segment_nodule"),
                       filter=formatter(),
                       output="{path[0]}/{basename[0]}-" + iso_size + "mm.nrrd")

    pipeline.transform(name="task_feature_extraction_seg",
                       task_func=tasks.feature_extraction,
                       input=output_from("task_image_resample_seg"),
                       filter=formatter(
                           experiment_set + r"/(?P<pid>[^/]+)/(?P<image_name>.*)_(?P<nid>\d*)-seg-label-" + iso_size + "mm.nrrd"),
                       output=output_path + "/{image_name[0]}_{nid[0]}-seg-" + iso_size + "mm.txt",
                       extras=["{subpath[0][1]}/CT/{image_name[0]}.nrrd"]) \

    pipeline.merge(name="task_feature_organization",
                   task_func=tasks.feature_organization,
                   input=output_from(["task_feature_extraction", "task_feature_extraction_seg"]),
                   output=feature_list_path)

    return pipeline


if __name__ == "__main__":
    freeze_support()

    pipeline_lungx_train = make_pipeline_lungx(experiment_set)

    parser = cmdline.get_argparse(description='TrainingSet lungx radiomics')

    options = parser.parse_args()

    # standard python logger which can be synchronised across concurrent
    # Ruffus tasks
    logger, logger_mutex = cmdline.setup_logging(
        __name__, options.log_file, options.verbose)

    #
    #   Run
    #
    cmdline.run(options)

