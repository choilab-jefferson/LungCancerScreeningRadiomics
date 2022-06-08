# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:55:23 2016

@author: choiw
"""

import pandas as pd
import numpy as np
import scipy.spatial.distance as spdist


class Nodules:
    nodules = pd.DataFrame()
    combined_nodules = pd.DataFrame()

    def __init__(self, filename=""):
        if len(filename) > 0:
            self.load_nodule_info(filename)

    def load_nodule_info(self, filename):
        self.nodules = pd.read_csv(filename)

    def find_same_nodules(self):
        nodules = self.nodules
        nodules.insert(nodules.columns.get_loc('nid'), column='nnid', value=0)
        nodules.insert(nodules.columns.get_loc('nid'), column='phys', value=0)

        distm = spdist.squareform(spdist.pdist(
            nodules[['Centroid_idx_1', 'Centroid_idx_2', 'Centroid_idx_3', ]]))
        for row in range(distm.shape[0]):
            th = np.average(
                nodules[['BoundingBox_idx_4', 'BoundingBox_idx_5']].iloc[row]) / 4
            count = 0
            for col in range(row, distm.shape[1]):
                if distm[row, col] < th and nodules['nnid'][col] == 0:  # find similar centroids
                    count = count + 1
                    nodules.at[col, 'nnid'] = row + 1000
                    nodules.at[col, 'phys'] = count
                    nodules.at[col, 'sid'] = count
                    print(row, col, distm[row, col], th, row + 1000)

        nodules.sort_values(['Volume', 'phys'], ascending=False, inplace=True)

        # re-numbering
        nnids = []
        nnnid = 1
        onnid = nodules['nnid'].copy()
        for nnid in onnid:
            if nnid not in nnids:
                nodules.loc[nodules['nnid'] == nnid, 'nnid'] = nnnid
                nnids.append(nnid)
                print(nnid, nnids)
                nnnid = nnnid + 1
        nodules.rename(columns={'nid': 'onid', 'nnid': 'nid'}, inplace=True)
        self.nodules = nodules.sort_values(['sid', 'nid'])

        return nodules

    def combine_same_nodules(self):
        nodules = self.nodules
        nodules.drop(['sid', 'onid'], axis=1, inplace=True)
        columns = nodules.columns
        dtypes = nodules.dtypes[columns]
        nids = np.unique(nodules['nid'])
        pid = nodules['pid'][0]
        new_nodules = []
        for nid in nids:
            selected = nodules['nid'] == nid

            #Boundary = nodules.loc[selected,['Lobulation','Spiculation']]
            #mean_Boundary = Boundary[Boundary>1].mean().round()
            # mean_Boundary[np.isnan(mean_Boundary)]=1
            #nodules.loc[selected,['Lobulation','Spiculation']]= (Boundary>1).astype(int)

            info_idx = range(columns.get_loc('Volume'),
                             columns.get_loc('BoundingBox_idx_6') + 1)
            ch_idx = range(columns.get_loc('Subtlety'),
                           columns.get_loc('Malignancy') + 1)
            bb_idx1 = range(columns.get_loc('BoundingBox_idx_1'),
                            columns.get_loc('BoundingBox_idx_3') + 1)
            bb_idx2 = range(columns.get_loc('BoundingBox_idx_4'),
                            columns.get_loc('BoundingBox_idx_6') + 1)
            phys = nodules.loc[selected, 'phys'].max()
            mean_info = nodules.loc[selected, nodules.columns[info_idx]].mean()
            print((np.asarray(bb_idx1) - info_idx[0]).tolist(), info_idx[
                  0], mean_info.index.get_loc('BoundingBox_idx_1'), nodules.loc[selected, nodules.columns[bb_idx1]].min())
            mean_info.iloc[(np.asarray(bb_idx1) - info_idx[0])
                           ] = np.floor(nodules.loc[selected, nodules.columns[bb_idx1]].min())
            mean_info.iloc[(np.asarray(bb_idx2) - info_idx[0])
                           ] = np.ceil(nodules.loc[selected, nodules.columns[bb_idx2]].max())
            mean_info[['MinIntensity', 'MaxIntensity']] = np.round(
                mean_info[['MinIntensity', 'MaxIntensity']])
            min_ch = nodules.loc[selected, nodules.columns[ch_idx]].min()
            max_ch = nodules.loc[selected, nodules.columns[ch_idx]].max()
            range_ch = max_ch - min_ch
            mode_ch = nodules.loc[selected, nodules.columns[ch_idx]].mode()
            mean_ch = nodules.loc[selected, nodules.columns[ch_idx]].astype('float').mean().round()

            nodule = {'pid': pid, 'nid': nid, 'phys': phys, **mean_info, **mean_ch}
            #print(nodules.iloc[selected,])
            for idx in ch_idx:
                idx1 = idx - columns.get_loc('Subtlety')
                if range_ch[idx1] > 2 and mode_ch.shape[0]>0 and not np.isnan(mode_ch.iloc[0,idx1]):
                    nodule[columns[idx]] = mode_ch.iloc[0,idx1]
                print(columns[idx], nodule[columns[idx]], range_ch[idx1], max_ch[idx1])

            if nodule['Texture'] == 2:
                nodule['Texture'] = 3
            if nodule['Texture'] == 4:
                nodule['Texture'] = 5

            new_nodules.append(nodule)
        new_nodules = pd.DataFrame(new_nodules, columns=columns)
        new_nodules = new_nodules.astype(dtypes, copy=False)
        self.combined_nodules = new_nodules
        return self.combined_nodules

    def save_combined_noduels(self, filename):
        self.combined_nodules.to_csv(filename, index=False)
