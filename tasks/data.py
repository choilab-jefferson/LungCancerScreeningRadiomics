import os

import os.path as osp
import numpy as np
import pandas as pd
import SimpleITK as sitk
import subprocess

from utils import nodule_info
from utils import LungRADS

import pylidc as pl
from multiprocessing import Manager

try:
    manager = Manager()
    mutex = manager.Lock()
except:
    pass

def dicom_to_nrrd_convert(input_dicom_path, output_file):
    dicom2nrrd_converter = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools DICOM-RT2NRRDConverter".split(" ")
    image_resample = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools python PythonTools/image_resample.py".split(" ")
    input_dicom_path = osp.abspath(input_dicom_path)
    output_file = osp.abspath(output_file)
    print('convert', input_dicom_path, output_file)
    for d in os.listdir(input_dicom_path):
        if d[0] != '.':
            series = d
    for d in os.listdir(os.path.join(input_dicom_path, series)):
        if d[0] != '.':
            study = d
    input_file = os.path.join(input_dicom_path, series, study)
    image_file = output_file.replace(".nrrd", "-1mm.nrrd")
    p = subprocess.Popen(dicom2nrrd_converter + [input_file, 'no', output_file.replace(".nrrd", "")])
    p.wait()
    print(output_file, image_file)
    p = subprocess.Popen(image_resample + [output_file, image_file, "1", "1", "1"])
    p.wait()


def load_scan_list(starting_data_path, selected=None):
    ppid = ''
    ct_paths = []
    print("Load LIDC dataset")
    for scan in pl.query(pl.Scan):
        pid = scan.patient_id
        if ppid == pid:
            pid = pid + "-1"
        ppid = pid
        if selected and pid not in selected:
            continue
        if len(scan.annotations) == 0:
            continue
        ct_paths.append(osp.join(starting_data_path, "CT", f"{pid}_CT.nrrd"))
        
    print(f"Loaded {len(ct_paths)} scans")
    return ct_paths


def load_patients_list(starting_data_path):
    # retrieve nrrd files and preparation of the input image sets
    CT_path = osp.join(starting_data_path, "CT")
    org_patients_list = [osp.join(CT_path, fn) for fn in next(os.walk(CT_path))[2]]
    #print(org_patients_list)

    for patient in org_patients_list:
        pid = osp.basename(patient)
        #print(pid)
        if pid.find('LIDC') < 0:
            continue
        
        if pid[0] == '.' or pid.find('mm') > 0 or pid.find('label') > 0 or pid.find('csv') > 0 or pid.find('txt') > 0 or pid.find('levelset') > 0:
            continue

        pt_nrrd_path = osp.join(CT_path, pid)
        yield pt_nrrd_path


def originate(output_file):
    dirname = osp.dirname(output_file)
    dirname = osp.dirname(dirname)
    pid = osp.basename(output_file).split("_")[0]
    dirname = osp.join(dirname, pid)

    input_nodule_info = osp.join(dirname, pid + '.csv')
    output_nodule_info = osp.join(dirname, pid + '_all.csv')
    nodules = nodule_info.Nodules(input_nodule_info)
    nodules.find_same_nodules()
    new_nodules = nodules.combine_same_nodules()
    lungRADS = LungRADS()
    lungRADS.nodules = new_nodules
    lungRADS.estimate()
    nodules.combined_nodules = lungRADS.nodules
    nodules.save_combined_noduels(output_nodule_info)


def originate_pylidc(output_file):
    def write_volume_to_nrrd(volume, origin, spacing, output_file):
        image = sitk.GetImageFromArray(volume.transpose([2, 0, 1])) # xyz -> zxy
        image.SetOrigin(origin)
        image.SetSpacing(spacing)
        sitk.WriteImage(image, output_file, True)
    

    pid = osp.basename(output_file).split("_")[0]
    dirname = osp.join(osp.dirname(osp.dirname(output_file)), pid)
    print(pid)

    with mutex:
        try:
            scan = pl.query(pl.Scan).filter(pl.Scan.patient_id==pid).first()
            annotations = scan.annotations
        except: # second
            scan = pl.query(pl.Scan).filter(pl.Scan.patient_id==pid.replace("-1","")).all()[1]
            annotations = scan.annotations

    # load CT DICOM images and convert them to nrrd volume
    images = scan.load_all_dicom_images(False)
    volume = np.stack(
            [
                x.pixel_array * x.RescaleSlope + x.RescaleIntercept
                for x in images
            ],
            axis=-1,
        ).astype(np.int16)
    origin = images[0].ImagePositionPatient
    spacing = [*images[0].PixelSpacing,
            images[1].ImagePositionPatient[2]-images[0].ImagePositionPatient[2]]
    write_volume_to_nrrd(volume, origin, spacing, output_file)

    # save nodule annotation into csv
    input_nodule_info = osp.join(dirname, f"{pid}.csv") # not from this
    output_nodule_info = osp.join(dirname,f"{pid}_all.csv")
    print(f"  {pid} has {len(annotations)} nodules.")

    new_nodules = []
    for i, ann in enumerate(annotations):
        with mutex:
            nid = ann.id
            bbox = ann.bbox()

        new_nodule = dict(
            pid=pid,
            sid=0,
            nid=nid,
            Volume=ann.volume,
            FilledVolume=ann.volume,
            SolidVolume=ann.volume,
            MeanIntensity=0,
            MinIntensity=0,
            MaxIntensity=0,
            Centroid_idx_1=ann.centroid[1], # matlab order
            Centroid_idx_2=ann.centroid[0],
            Centroid_idx_3=ann.centroid[2],
            WeightedCentroid_idx_1=0,
            WeightedCentroid_idx_2=0,
            WeightedCentroid_idx_3=0,
            BoundingBox_idx_1=bbox[1].start,
            BoundingBox_idx_2=bbox[0].start,
            BoundingBox_idx_3=volume.shape[2]-bbox[2].stop,
            BoundingBox_idx_4=bbox[1].stop-bbox[1].start,
            BoundingBox_idx_5=bbox[0].stop-bbox[0].start,
            BoundingBox_idx_6=bbox[2].stop-bbox[2].start,
            Subtlety=ann.subtlety,
            InternalStructure=ann.internalStructure,
            Calcification=ann.calcification,
            Sphericity=ann.sphericity,
            Margin=ann.margin,
            Lobulation=ann.lobulation,
            Spiculation=ann.spiculation,
            Texture=ann.texture,
            Malignancy=ann.malignancy,
        )
        new_nodules.append(new_nodule)

    os.makedirs(dirname, exist_ok=True)
    nodules = nodule_info.Nodules()
    nodules.nodules = pd.DataFrame(new_nodules)
    nodules.find_same_nodules()
    nodules_df = nodules.nodules.copy()
    nodules_df.to_csv(input_nodule_info, index=False) #
    nodules.combine_same_nodules()
    lungRADS = LungRADS()
    lungRADS.nodules = nodules.combined_nodules
    lungRADS.estimate()
    nodules.combined_nodules = lungRADS.nodules
    nodules.save_combined_noduels(output_nodule_info)

    # extract nodule masks
    nodule_label = np.zeros(volume.shape, np.ubyte)
    nodule_label_file = osp.join(dirname, f"{pid}_CT-label.nrrd")
    for si in range(4):
        si_nodules_df = nodules_df.loc[nodules_df.sid==si+1,]
        sid = 'Phy' + str(si + 1)
        si_nodule_label_file = osp.join(dirname, f"{pid}_CT_{sid}-label.nrrd")
        si_factor = 2**si
        si_nodule_label = np.zeros(volume.shape, np.ubyte)
        for ni, nodule in si_nodules_df.iterrows():
            with mutex:
                ann = pl.query(pl.Annotation).filter(pl.Annotation.id == nodule.onid).first()
                mask = ann.boolean_mask()
                bbox = ann.bbox()
            si_nodule_label[bbox] += (mask).astype(np.ubyte)
        write_volume_to_nrrd(si_nodule_label, origin, spacing, si_nodule_label_file)
        nodule_label += si_nodule_label*si_factor
    write_volume_to_nrrd(nodule_label, origin, spacing, nodule_label_file)
    
    print(output_file)


def check_nodules(input_file, output_file):
    print(input_file)
    print(output_file)
    dirname = osp.dirname(input_file)
    dirname = osp.dirname(dirname)
    pid = osp.basename(input_file).split("_")[0]
    dirname = osp.join(dirname, pid)

    input_nodule_info = osp.join(dirname, pid + '_all.csv')
    nodules = nodule_info.Nodules(input_nodule_info).nodules

    for nidx, nodule in nodules.iterrows():
        #if nidx != 0:
        #    continue
        with open(osp.join(dirname, pid + '_CT_' + str(nidx + 1) + '.nodule'), "w"):
            pass

