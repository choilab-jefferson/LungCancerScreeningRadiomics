#!/usr/bin/env python3
from tciaexplorer import TciaExplorer
import os
import json
import requests
import sys
from zipfile import ZipFile
from tqdm import tqdm
import pandas as pd


def download_all_images(collection, location):
    tcia = TciaExplorer()
    
    print ("Fetching patients for collection: "+ collection)
    ############Fetch Patients#############
    try:
        patients = json.loads(tcia.get_patient(collection=collection).text)
        patients = pd.DataFrame(patients).PatientID
    except requests.exceptions.RequestException as e:
        print(e)
        return -1

    for patient in tqdm(patients):
        ###################Fetch study###############
        print ("Fetching study for patientID: "+patient)
        try:
            patientStudy  = json.loads(tcia.get_patient_study(collection=collection, patientID=patient).text)
            print(patientStudy)
        except requests.exceptions.RequestException as e:
            print(e)
            return -1

        for study in patientStudy:
            ##########Fetch series#############
            print("Fetching series for the studyInstanceUID: "+ study['StudyInstanceUID'])
            try:
                patientSeries = json.loads(tcia.get_series(studyInstanceUID=study['StudyInstanceUID']).text)
            except requests.exceptions.RequestException as e:
                print(e)
                return -1

            for series in patientSeries:
                print("Fetching images for the seriesInstanceUID: "+series['SeriesInstanceUID'])
                ##################Fetch images#############
                try:
                    images = (tcia.get_image(seriesInstanceUID=series['SeriesInstanceUID']))
                except requests.exceptions.RequestException as e:
                    print(e)
                    return -1

                #print(images)
                print("Writing image zipfile for the seriesInstanceUID: "+series['SeriesInstanceUID'])
                series_path = os.path.join(location, patient,study['StudyInstanceUID'], series['SeriesInstanceUID'])

                os.makedirs(series_path, exist_ok=True)

                fileName = os.path.join(series_path,patient+".zip")
                f = open(fileName, "wb")
                f.write(images.content)
                f.close()

                print("Unzipping images for the seriesInstanceUID: "+series['SeriesInstanceUID'])
                with ZipFile(fileName, 'r') as myzip:
                    myzip.extractall(series_path)
                os.unlink(fileName)

    return 0

if __name__ == "__main__":
    download_all_images('SPIE-AAPM Lung CT Challenge', 'DATA/LUNGx')
    download_all_images('LIDC-IDRI', 'DATA/LIDC-IDRI')
