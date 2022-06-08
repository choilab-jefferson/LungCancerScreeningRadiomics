import os
import sys
import subprocess

import os.path as osp
import numpy as np
import SimpleITK as sitk

def _image_resample(input_image, ref_image, is_label=False):
    resample_filter = sitk.ResampleImageFilter()

    input_spacing = input_image.GetSpacing()
    input_direction = input_image.GetDirection()
    input_origin = input_image.GetOrigin()
    input_size = input_image.GetSize()

    output_spacing = ref_image.GetSpacing()
    output_direction = ref_image.GetDirection()

    flip = ((input_direction[0] > 0) != (output_direction[0] > 0),
            (input_direction[4] > 0) != (output_direction[4] > 0),
            (input_direction[8] > 0) != (output_direction[8] > 0))

    input_image = sitk.Flip(input_image, flip)
    input_image.SetDirection(output_direction)

    input_direction = input_image.GetDirection()
    input_origin = input_image.GetOrigin()

    #print(flip, (input_direction[0] > 0, output_direction[0] > 0,
    #             input_direction[4] > 0, output_direction[4] > 0,
    #             input_direction[8] > 0, output_direction[8] > 0))

    origin = ref_image.GetOrigin()
    dist = (np.asarray(input_origin) - origin) / np.asarray(input_spacing)
    close_interp = np.floor(
        dist * np.asarray(input_spacing) / np.asarray(output_spacing))
    output_origin = ref_image.GetOrigin() + close_interp * \
        np.asarray(output_spacing)

    output_size = np.ceil(np.asarray(
        input_size) * np.asarray(input_spacing) / np.asarray(output_spacing)).astype(int)

    resample_filter.SetOutputSpacing(output_spacing)
    resample_filter.SetOutputOrigin(output_origin)
    resample_filter.SetSize(output_size.tolist())
    resample_filter.SetOutputDirection(output_direction)
    if is_label:
        resample_filter.SetInterpolator(sitk.sitkNearestNeighbor)
        output_image = resample_filter.Execute(input_image)
    else:
        resample_filter.SetInterpolator(sitk.sitkLinear)
        output_image = resample_filter.Execute(input_image)

    resampled_image = resample_filter.Execute(input_image)

    return resampled_image


def _compare_masks(image1, image2):
    overlap_measures = sitk.LabelOverlapMeasuresImageFilter()
    overlap_measures.Execute(image1, image2)
    print(overlap_measures)

    return overlap_measures


def _merge_gc_cip(input_files, output_file):
    gc_mask_file = input_files[0]
    cip_mask_file = input_files[1]
    gc_mask_image = sitk.ReadImage(gc_mask_file)
    cip_mask_image = sitk.ReadImage(cip_mask_file)

    # resample cip to match gc
    cip_mask_image = _image_resample(cip_mask_image, gc_mask_image, True) # to make it in the same space
    cip_mask_image = sitk.ConstantPad(cip_mask_image, (100, 100, 100), (100, 100, 100))
    cip_mask_image = sitk.RegionOfInterest(cip_mask_image, gc_mask_image.GetSize(),cip_mask_image.TransformPhysicalPointToIndex(gc_mask_image.GetOrigin()))
    cip_mask_image.SetOrigin(gc_mask_image.GetOrigin())
    cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    # cip_mask_image = sitk.BinaryClosingByReconstruction(cip_mask_image, 1, sitk.sitkBall, 1, True)
    # cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    cip_mask_image = sitk.Cast(cip_mask_image, sitk.sitkInt32)

    sitk.WriteImage(cip_mask_image, output_file.replace('seg', 'cip-iso'), True)

    ##    
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

    if G_ > C_:
        bigger_mask_image,smaller_mask_image = gc_mask_image, cip_mask_image
    else:
        bigger_mask_image,smaller_mask_image = cip_mask_image, gc_mask_image
    if (G_<=0 or C_<=0) and I/U<0.5:
        smaller_mask_image = bigger_mask_image

    
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
    sitk.WriteImage(attched_image, output_file.replace('seg', 'attached'), True)


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
    sitk.WriteImage(out_mask_image, output_file.replace('seg', 'seg-iso'), True)

    
def _resample_segmentations(input_files, output_file):
    gc_mask_file = input_files[0]
    cip_mask_file = input_files[1]
    gc_mask_image = sitk.ReadImage(gc_mask_file)
    cip_mask_image = sitk.ReadImage(cip_mask_file)

    ## compare with physician's contour
    aniso_gc_mask_image = _image_resample(gc_mask_image, mask_image, True)
    aniso_cip_mask_image = _image_resample(cip_mask_image, mask_image, True)

    mask_image = sitk.ConstantPad(mask_image, (100, 100, 100), (100, 100, 100))
    mask_image = sitk.RegionOfInterest(mask_image, aniso_gc_mask_image.GetSize(),
        mask_image.TransformPhysicalPointToIndex(aniso_gc_mask_image.GetOrigin()))
    mask_image.SetOrigin(aniso_gc_mask_image.GetOrigin())

    aniso_out_mask_image = _image_resample(out_mask_image, mask_image, True)
    aniso_out_mask_image.SetOrigin(mask_image.GetOrigin())

    sitk.WriteImage(aniso_gc_mask_image, output_file.replace('seg', 'gc'), True)
    sitk.WriteImage(aniso_cip_mask_image, output_file.replace('seg', 'cip'), True)
    sitk.WriteImage(aniso_out_mask_image, output_file, True)


def _resample_segmentations(input_files, output_file):
    gc_mask_file = input_files[0]
    cip_mask_file = input_files[1]
    cip_mask_file = input_files[3]
    gc_mask_image = sitk.ReadImage(gc_mask_file)
    cip_mask_image = sitk.ReadImage(cip_mask_file)

    ##    
    union_image = sitk.Or(cip_mask_image, gc_mask_image)
    intersection_image = sitk.And(cip_mask_image, gc_mask_image)

    ## select bigger and smaller segmentation from growcut and cip
    U = (sitk.GetArrayFromImage(union_image)>0).sum()
    I = (sitk.GetArrayFromImage(intersection_image)>0).sum()
    G = (sitk.GetArrayFromImage(gc_mask_image)>0).sum()
    C = (sitk.GetArrayFromImage(cip_mask_image)>0).sum()
    G_ = G-I
    C_ = C-I

    gc_compare = _compare_masks(mask_image, aniso_gc_mask_image)
    cip_compare = _compare_masks(mask_image, aniso_cip_mask_image)
    seg_compare = _compare_masks(mask_image, aniso_out_mask_image)

    fsegmesure = open(mask_file.replace('-label.nrrd', '.txt'),'w')
    print(U,I,G,C,I/U,G/U,C/U,G_,C_, file=fsegmesure)
    print("Growcut", file=fsegmesure)
    print(gc_compare, file=fsegmesure)
    print("CIP", file=fsegmesure)
    print(cip_compare, file=fsegmesure)
    print("Segmentation", file=fsegmesure)
    print(seg_compare, file=fsegmesure)
    fsegmesure.close()

def _growcut_segmentation_with_mask(image_file, gc_mask_file, in_mask_file):
    input_image = sitk.ReadImage(image_file)
    mask_image = sitk.ReadImage(in_mask_file)

    filter = sitk.LabelShapeStatisticsImageFilter()
    filter.ComputeFeretDiameterOn()
    filter.Execute(mask_image)
    center = filter.GetCentroid(1)
    radius = filter.GetFeretDiameter(1)/2
    #print("center", center)
    #print("radius", radius)
    center_idx = input_image.TransformPhysicalPointToIndex(center)

    ## grow cut segmenation
    growcut_segmentation(image_file, gc_mask_file, center_idx, radius)
    
    # post processing
    gc_mask_image = sitk.ReadImage(gc_mask_file)
    gc_mask_image = sitk.Cast(sitk.BinaryThreshold(gc_mask_image, 1, 1000), sitk.sitkInt32)
    gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)
    # gc_mask_image = sitk.BinaryClosingByReconstruction(gc_mask_image, 1, sitk.sitkBall, 1, True)
    # gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)

    sitk.WriteImage(gc_mask_image, gc_mask_file, True)


def _cip_segmentation_with_mask(input_files, output_file, in_mask_file):
    image_file = input_files[0]
    image_crop_file = input_files[0]
    
    mask_image = sitk.ReadImage(in_mask_file)
    input_image_crop = sitk.ReadImage(image_crop_file)

    filter = sitk.LabelShapeStatisticsImageFilter()
    filter.ComputeFeretDiameterOn()
    filter.Execute(mask_image)
    center = filter.GetCentroid(1)
    radius = filter.GetFeretDiameter(1)/2
    
    ## chest image platform (CIP) segmenation
    partSolid = '--echo'
    # if nodule.loc['Texture'] != 5:
    #    partSolid = '--partSolid'

    # seed generation
    seeds = list()
    seeds.append('--seeds')
    seeds.append('%f,%f,%f' % center)
    try:
        np_mask_image = sitk.GetArrayFromImage(mask_image) > 0
        np_input_image = sitk.GetArrayFromImage(input_image_crop)

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

    levelset_file = output_file.replace("cip","levelset")
    _cip_segmentation(image_file, levelset_file, radius, seeds, partSolid)
    
    try:
        levelset_image = sitk.ReadImage(levelset_file)
    except:
        raise('Failed CIP')

    ## post processing
    cip_mask_image = sitk.BinaryThreshold(levelset_image, 0)
    sitk.WriteImage(cip_mask_image, output_file, True)


def _growcut_segmentation(image_file, gc_mask_file, center_idx, radius):
    nodule_segmentation = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools NoduleSegmentation".split(" ")
    image_file = osp.abspath(image_file)
    gc_mask_file = osp.abspath(gc_mask_file)
    p = subprocess.Popen(nodule_segmentation + [image_file, str(center_idx[0]),
                        str(center_idx[1]), str(-center_idx[2]), str(radius*2+5), str(radius), gc_mask_file])
    p.wait()


def _cip_segmentation(image_file, levelset_file, radius, seeds, partSolid):
    nodule_segmentation_cip = f"docker run --volume {os.getcwd()}:{os.getcwd()} acilbwh/chestimagingplatform GenerateLesionSegmentation".split(" ")
    image_file = osp.abspath(image_file)
    levelset_file = osp.abspath(levelset_file)
    p = subprocess.Popen(
        nodule_segmentation_cip + ['-i', image_file, '-o', levelset_file, '--maximumRadius', str(radius+5), partSolid, '--echo'] + seeds)
    p.wait()


def segment_nodule(input_files, output_files, image_file):
    in_mask_file, staple_file = input_files[0], input_files[1]
    mask_file, levelset_file = output_files[0], output_files[1]

    if False:
        dirname = osp.dirname(in_mask_file)
        names = osp.basename(in_mask_file).split("_")
        basename = names[0]
        input_nodule_info = osp.join(dirname, basename + '_all.csv')
        nidx = int(names[2].replace('-all-label.nrrd', '')) - 1

        nodules = nodule_info.Nodules(input_nodule_info).nodules
        nodule = nodules.loc[nidx, ]


    input_image = sitk.ReadImage(image_file)
    staple_image = sitk.ReadImage(staple_file)
    mask_image = sitk.ReadImage(in_mask_file)

    input_image_crop = sitk.RegionOfInterest(input_image, mask_image.GetSize(), input_image.TransformPhysicalPointToIndex(mask_image.GetOrigin()))
    sitk.WriteImage(input_image_crop, in_mask_file.replace('-label', ''), True)
    
    np_staple_image = sitk.GetArrayFromImage(staple_image)
    np_mask_image = sitk.GetArrayFromImage(mask_image) > 0
    np_input_image = sitk.GetArrayFromImage(input_image_crop)

    filter = sitk.LabelShapeStatisticsImageFilter()
    filter.ComputeFeretDiameterOn()
    filter.Execute(mask_image)
    center = filter.GetCentroid(1)
    radius = filter.GetFeretDiameter(1)/2
    print("center", center)
    print("radius", radius)

    
    print(mask_image.GetSize())
    print(input_image_crop.GetSize())
    
    ## grow cut segmenation
    gc_mask_file = mask_file.replace('seg', 'growcut')
    center_idx = input_image.TransformPhysicalPointToIndex(center)
    print(center_idx)
    _growcut_segmentation(image_file, gc_mask_file, center_idx, radius)
    gc_mask_image = sitk.ReadImage(gc_mask_file)
    gc_mask_image = sitk.Cast(sitk.BinaryThreshold(gc_mask_image, 1, 1000), sitk.sitkInt32)
    
    # post processing
    gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)
    # gc_mask_image = sitk.BinaryClosingByReconstruction(
    #     gc_mask_image, 1, sitk.sitkBall, 1, True)
    # gc_mask_image = sitk.VotingBinaryIterativeHoleFilling(gc_mask_image)

    sitk.WriteImage(gc_mask_image, mask_file.replace('seg', 'gc-iso'), True)
    
    
    ## make physician's contour to isotropic voxel
    mask_image = sitk.Cast(sitk.BinaryThreshold(mask_image, 1, 1000), sitk.sitkInt32)
    iso_mask_image = _image_resample(mask_image, gc_mask_image, True)
    iso_mask_image = sitk.ConstantPad(iso_mask_image, (100, 100, 100), (100, 100, 100))
    iso_mask_image = sitk.RegionOfInterest(iso_mask_image, gc_mask_image.GetSize(),
                       iso_mask_image.TransformPhysicalPointToIndex(gc_mask_image.GetOrigin()))
    iso_mask_image.SetOrigin(iso_mask_image.GetOrigin())
    sitk.WriteImage(iso_mask_image, mask_file.replace('seg', 'all-iso'), True)

    ## chest image platform (CIP) segmenation
    partSolid = '--echo'
    # if nodule.loc['Texture'] != 5:
    #    partSolid = '--partSolid'


    # seed generation
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

    _cip_segmentation(image_file, levelset_file, radius, seeds, partSolid)

    try:
        levelset_image = sitk.ReadImage(levelset_file)
    except:
        print('Failed CIP')
        levelset_image = iso_mask_image

    if (sitk.GetArrayFromImage(gc_mask_image)>0).sum() < 9:
        print('Failed GC')
        gc_mask_image = iso_mask_image
    

    sys.stdout.flush()


    ## post processing
    cip_mask_image = sitk.BinaryThreshold(levelset_image, 0)
    cip_mask_image = _image_resample(cip_mask_image, gc_mask_image, True) # to make it in the same space
    cip_mask_image = sitk.ConstantPad(cip_mask_image, (100, 100, 100), (100, 100, 100))
    cip_mask_image = sitk.RegionOfInterest(cip_mask_image, gc_mask_image.GetSize(),
                                          cip_mask_image.TransformPhysicalPointToIndex(gc_mask_image.GetOrigin()))
    cip_mask_image.SetOrigin(gc_mask_image.GetOrigin())
    cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    # cip_mask_image = sitk.BinaryClosingByReconstruction(
    #     cip_mask_image, 1, sitk.sitkBall, 1, True)
    # cip_mask_image = sitk.VotingBinaryIterativeHoleFilling(cip_mask_image)
    cip_mask_image = sitk.Cast(cip_mask_image, sitk.sitkInt32)


    sitk.WriteImage(cip_mask_image, mask_file.replace('seg', 'cip-iso'), True)


    ##    
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

    if G_ > C_:
        bigger_mask_image,smaller_mask_image = gc_mask_image, cip_mask_image
    else:
        bigger_mask_image,smaller_mask_image = cip_mask_image, gc_mask_image
    if (G_<=0 or C_<=0) and I/U<0.5:
        smaller_mask_image = bigger_mask_image

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

    
    ## compare with physician's contour
    aniso_gc_mask_image = _image_resample(gc_mask_image, mask_image, True)
    aniso_cip_mask_image = _image_resample(cip_mask_image, mask_image, True)

    mask_image = sitk.ConstantPad(mask_image, (100, 100, 100), (100, 100, 100))
    mask_image = sitk.RegionOfInterest(mask_image, aniso_gc_mask_image.GetSize(),
        mask_image.TransformPhysicalPointToIndex(aniso_gc_mask_image.GetOrigin()))
    mask_image.SetOrigin(aniso_gc_mask_image.GetOrigin())

    aniso_out_mask_image = _image_resample(out_mask_image, mask_image, True)
    aniso_out_mask_image.SetOrigin(mask_image.GetOrigin())

    sitk.WriteImage(aniso_gc_mask_image, mask_file.replace('seg', 'gc'), True)
    sitk.WriteImage(aniso_cip_mask_image, mask_file.replace('seg', 'cip'), True)
    sitk.WriteImage(aniso_out_mask_image, mask_file, True)

    gc_compare = _compare_masks(mask_image, aniso_gc_mask_image)
    cip_compare = _compare_masks(mask_image, aniso_cip_mask_image)
    seg_compare = _compare_masks(mask_image, aniso_out_mask_image)


    fsegmesure = open(mask_file.replace('-label.nrrd', '.txt'),'w')
    print(U,I,G,C,I/U,G/U,C/U,G_,C_, file=fsegmesure)
    print("Growcut", file=fsegmesure)
    print(gc_compare, file=fsegmesure)
    print("CIP", file=fsegmesure)
    print(cip_compare, file=fsegmesure)
    print("Segmentation", file=fsegmesure)
    print(seg_compare, file=fsegmesure)
    fsegmesure.close()