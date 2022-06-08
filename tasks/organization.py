import os

import pandas as pd

import sys
from utils import organize_features

def feature_organization(input_files, output_file):
    organize_features.organize(
        input_files, output_file, organize_features.parameter_list, organize_features.feature_list)


def staple_organization(input_files, output_file):
    staple_result_list = []
    for input_file in input_files:
        staple_result_file = input_file[2]
        staple_result = pd.read_csv(staple_result_file, sep='\t', header=0)
        l = [i.split('_') for i in staple_result['File ']]
        staple_result.drop(['File '], axis=1, inplace=True)
        info = pd.DataFrame(l, columns=['PID', 'Image', 'Tags'])
        info.iloc[len(l) - 1, 0] = info.iloc[0, 0]
        info.iloc[len(l) - 1, 1] = info.iloc[0, 1]
        info.iloc[len(l) - 1, 2] = info.iloc[0, 2].split('-')[0] + '-All'
        staple_result = pd.concat([info, staple_result], axis=1)
        staple_result_list.append(staple_result)
        #print(staple_result)
    #print(staple_result_list)
    staple_result_list = pd.concat(staple_result_list)
    staple_result_list.to_csv(output_file)


def segmentation_organization(input_files, output_file):
    def find(inputlist, value):
        try:
            idx = inputlist.index(value)
        except ValueError:
            idx = -1

        return idx


    def measurement_parsing(filename, in_measurement_list):
        len_m = len(in_measurement_list)
        values = [['']*len_m for i in range(3)]
        seg_idx = [0] * len_m
        with open(filename, 'r') as f:
            for line in f:
                strs = line.split(': ')
                if len(strs) == 2:
                    feature_name = strs[0].strip()
                    value_str = strs[1].strip()
                else:
                    continue

                value = value_str

                fidx = find(in_measurement_list, feature_name)
                if fidx > -1:
                    values[seg_idx[fidx]][fidx] = value
                    seg_idx[fidx] += 1
                    print(seg_idx[fidx], feature_name + "=" + value)

        return values

    segmentations = ['GrowCut', 'CIP', 'Final']
    measurements = ['DiceCoefficient','JaccardCoefficient','FalseNegativeError','FalsePositiveError','MeanOverlap','UnionOverlap','VolumeSimilarity']
    with open(output_file, 'w') as fw:
        column_names = ','.join(['No.', 'PID', 'Tags'] + measurements)
        fw.write(column_names + '\n')

        print(input_files)

        no = 1
        for f in input_files:
            f = f[0].replace('-label.nrrd','.txt')

            print(f)
            path = os.path.dirname(f)
            strs = os.path.basename(f).split("_")
            pid = strs[0]
            image = strs[1]

            values_matrix = measurement_parsing(f, measurements)
            print(values_matrix)

            for v in range(3):
                values_list = values_matrix[v]
                fw.write(str(no) + ',')
                fw.write(pid + ',')
                fw.write(segmentations[v] + ',')

                values = ','.join(values_list)
                fw.write(values+'\n')
                no += 1
