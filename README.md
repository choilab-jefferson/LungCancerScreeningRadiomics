# Lung Cancer Screening Radiomics

## Datasets
LIDC-IDRI https://wiki.cancerimagingarchive.net/display/Public/LIDC-IDRI  
LUNGx https://wiki.cancerimagingarchive.net/display/Public/LUNGx+SPIE-AAPM-NCI+Lung+Nodule+Classification+Challenge  

Download all the required datasets from TCIA using the following command. 
```bash
$ python3 download_datasets.py
```
This will store the datasets into `DATA/LIDC-IDRI` and `DATA/LUNGx`.
This is a single threaded script. You better use TCIA downloder for the best speed in the links above.


## Usage
This radiomics pipeline uses docker to run radiomics-tools and other external tools.

### 1. Preprocessing and feature extractions

#### * LIDC

1. DICOM conversion (matlab)
   * input - DATA/LIDC-IDRI (download from TCIA)
   * output - DATA/nodule-lidc   
   ```matlab
   >> cd matlab
   >> run main_image_to_nrrd.m
   ```

2. Radiomics feature extraction (python)
   * input - DATA/LIDC-IDRI
   * output - DATA/nodule-lidc, output/nodule-lidc, output/feature-list_nodule-lidc.csv, output/STAPLE_nodule-lidc.csv  
   ```bash
   $ python3 run_lidc.py
   ```


Note: `$ python3 run_lidc_pylidc.py` is also available to extract radiomics features without matlab, but it cannot generate full nodule information list.


#### * LUNGx
   * input - DATA/LUNGx  
   * output - DATA/nodule-lungx, output/nodule-lungx, output/feature_list_lungx.csv
   ```bash
   $ python3 run_lungx.py
   ```


### 2. Spiculation feature extraction

1. Mesh model generation from manual or semi-automatic segmentation  (matlab)
   * input - DATA/nodule-\*  
   * output - DATA/nodule-\*/objs  
   ```matlab
   >> cd matlab
   >> run convert_for_sphere_param.m % for all LIDC data
   >> run convert_for_sphere_param_lungx.m  % for LUNGx data
   ```

2. Spherical parameterization (matlab)
   * input - DATA/nodule-\*/objs  
   * output - DATA/nodule-\*/spherical_objs
   ```matlab
   >> cd matlab
   >> run main_sphere_param.m % for all LIDC data
   >> run main_sphere_param_lungx.m  % for LUNGx data
   ```


3. Spiculation feature extraction (matlab)
   * input - DATA/nodule-\*/spherical_objs\*
   * output - output/nodule-\*/spiculation\*, output/nodule-\*/parameters
    ```matlab
    >> cd matlab
    >> run main_spiculation.m % for all LIDC data
    >> run main_spiculation_lungx.m  % for LUNGx data
    ```

### 3. Malignancy prediction model

1. AutoML (python)
  It is possible to build and validate simple radiomics models using AutoML.
  In Jupyter, open automl.ipynb and run all cells. You may also use a python script to test it, as seen below.  
    
   * input - meatadata, output/feature-list_nodule-lidc.csv, output/feature-list_nodule-lungx.csv
   * output - output/nodule-lidc/spherical_objs
  
    ```bash
    $ python automl.py
    ```

## Publications
1. Wookjin Choi, Jung Hun Oh, Sadegh Riyahi, Chia-Ju Liu, Feng Jiang, Wengen Chen, Charles White, Andreas Rimner, James G. Mechalakos, Joseph O. Deasy, and Wei Lu, “Radiomics analysis of pulmonary nodules in low-dose CT for early detection of lung cancer”, Medical Physics, Vol. 45, No. 4, pp. 1537-1549, April 2018. https://doi.org/10.1002/mp.12820
2. Wookjin Choi, Saad Nadeem, Sadegh Riyahi, Joseph O. Deasy, Allen Tannenbaum, Wei Lu, “Reproducible and Interpretable Spiculation Quantification for Lung Cancer Screening.” Computer methods and programs in biomedicine. 200 (2021). https://doi.org/10.1016/j.cmpb.2020.105839