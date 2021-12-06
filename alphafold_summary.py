#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 12:31:02 2021

@author: andrea
"""

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data_frame = pd.read_csv('summary.csv')

data_frame_top = data_frame[data_frame['iptm']>0.670]
data_frame_top = data_frame_top.drop_duplicates(subset = ["run"])

plt.clf()
fig, axs = plt.subplots(ncols=6,
                        figsize = (16, 4))

sns.histplot(data_frame['ptm'], bins=20, ax=axs[0])
axs[0].title.set_text('pTM score')


sns.histplot(data_frame['iptm'], bins=20, ax=axs[1])
axs[1].title.set_text('ipTM score')

sns.histplot(data_frame['plddt'], bins=20, ax=axs[2])
axs[2].title.set_text('plddt')

sns.scatterplot(data_frame['ptm'], data_frame['plddt'], ax=axs[3])
axs[3].title.set_text('plddt vs pTM')


sns.scatterplot(data_frame['iptm'], data_frame['plddt'], ax=axs[4])
axs[4].title.set_text('plddt vs ipTM')

sns.scatterplot(data_frame_top['iptm'], data_frame_top['plddt'], ax=axs[5])
# add annotations one by one with a loop
axs[5].title.set_text('best hits by iptm')
# add annotations one by one with a loop
for line in range(0,data_frame_top.shape[0]):
     axs[5].text(data_frame_top.iptm.iloc[line]+0.001, 
              data_frame_top.plddt.iloc[line],
              data_frame_top.run.iloc[line], 
              horizontalalignment='left',
              size='small', 
              color='black')
     

plt.suptitle('Alphfold run summary')

plt.savefig('summary.png')


sorted_data_frame = data_frame.sort_values(by=['iptm'])
sorted_data_frame.to_csv('summary_sorted.csv')
