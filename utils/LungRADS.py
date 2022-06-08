# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 11:20:51 2016

@author: choiw
"""

import math
import pandas as pd

class LungRADS():
    nodules = pd.DataFrame()

    def __init__(self, filename = ""):
        if len(filename) > 0:
            self.load_nodule_info(filename)

    def load_nodule_info(self, filename):
        self.nodules = pd.read_csv(filename)


    def estimate(self):
        # now only support baseline
        # TODO: follow up screening

        nodules = self.nodules
        new_nodules = []
        for _, nodule in nodules.iterrows():
            nodule.D = math.pow(nodule.Volume*3/4/math.pi,1./3)*2

            nodule["LungRADS"] = '0'
            print("Diameter :", nodule.D)
            if nodule.Calcification != 4 and nodule.Calcification != 6:
                nodule.LungRADS = '1'
            else:
                if nodule.Texture == 1:
                    print('GGO')
                    if nodule.D < 20:
                        nodule.LungRADS = '2'
                    else:
                        nodule.LungRADS = '3'
                elif nodule.Texture == 3:
                    print('Part solid')
                    nodule.SD = math.pow(nodule.SolidVolume*3/4/math.pi,1./3)*2
                    if nodule.D < 6:
                        nodule.LungRADS = '2'
                    elif nodule.SD < 6:
                        nodule.LungRADS = '3'
                    elif nodule.SD < 8:
                        nodule.LungRADS = '4A'
                    elif nodule.SD >= 8:
                        nodule.LungRADS = '4B'
                elif nodule.Texture == 5:
                    print('Solid')
                    if nodule.D < 6:
                        nodule.LungRADS = '2'
                    elif nodule.D < 8:
                        nodule.LungRADS = '3'
                    elif nodule.D < 15:
                        nodule.LungRADS = '4A'
                    elif nodule.D >= 15:
                        nodule.LungRADS = '4B'
            if nodule.LungRADS == '3' or nodule.LungRADS == '4':
                if nodule.Spiculation >= 2:
                    nodule.LungRADS = '4X'

            print("LungRADS Category:", nodule.LungRADS)
            print("LIDC Malignancy:", nodule.Malignancy)
            new_nodules.append(nodule.to_dict())
        new_nodules = pd.DataFrame(new_nodules, columns=nodule.index)
        print(new_nodules)
        self.nodules = new_nodules
