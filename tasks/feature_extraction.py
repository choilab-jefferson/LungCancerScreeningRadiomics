import os
import subprocess
import os.path as osp


feature_extractor = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools FeatureExtraction".split(" ")
iso_size = "1"

def feature_extraction(mask_file, output_file, image_file):
    image_resample = f"docker run --volume {os.getcwd()}:{os.getcwd()} wookjinchoi/radiomics-tools python PythonTools/image_resample.py".split(" ")
    if type(mask_file) == list:
        mask_file = mask_file[0]
    mask_file_resize = mask_file.replace("-label.nrrd", "-" + iso_size + "mm-label.nrrd")
    print(mask_file, mask_file_resize)
    image_file_resize = image_file.replace(".nrrd", "-" + iso_size + "mm.nrrd")
    print(image_file, image_file_resize)
    p = subprocess.Popen(image_resample + [mask_file, mask_file_resize, image_file_resize, "1"])
    p.wait()
    p = subprocess.Popen(feature_extractor + [image_file_resize, mask_file_resize, osp.abspath(output_file), '1', '2igscr', '64'])
    p.wait()


def feature_extraction_phy(mask_file, output_file, image_file):
    for mi in range(4):
        if osp.getsize(mask_file[mi]) > 0:
            feature_extraction(
                [mask_file[mi]], output_file[mi], image_file)
        else:
            with open(output_file[mi], "w"):
                pass