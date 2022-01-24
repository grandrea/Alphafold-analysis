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

data_frame_top = data_frame[data_frame['iptm']>0.75]
data_frame_top = data_frame_top.sort_values(by=['iptm'], ascending=False)
data_frame_top = data_frame_top.drop_duplicates(subset = ["run"])

plt.clf()
fig, axs = plt.subplots(ncols=3,
                        figsize = (16, 3))

sns.histplot(data_frame['ptm'], bins=20, ax=axs[0])
axs[0].title.set_text('pTM score')


sns.histplot(data_frame['iptm'], bins=20, ax=axs[1])
axs[1].title.set_text('ipTM score')

sns.histplot(data_frame['plddt'], bins=20, ax=axs[2])
axs[2].title.set_text('plddt')


plt.suptitle('Alphfold run summary')
plt.savefig('summary_histograms.png', dpi=300)

plt.clf()
fig, axs = plt.subplots(ncols=3,
                        figsize = (16, 3))
sns.scatterplot(data_frame['ptm'], data_frame['plddt'], ax=axs[0])
axs[0].title.set_text('plddt vs pTM')

sns.scatterplot(data_frame['iptm'], data_frame['plddt'], ax=axs[1])
axs[1].title.set_text('plddt vs ipTM')

sns.scatterplot(data_frame['ptm'], data_frame['iptm'], ax=axs[2])
axs[2].title.set_text('ipTM vs pTM')


plt.suptitle('Alphfold run summary')
plt.savefig('summary_scatter.png', dpi=300)

plt.clf()
fig, axs = plt.subplots(ncols=2,
                        figsize = (12, 2))

sns.scatterplot(data_frame_top['iptm'], data_frame_top['plddt'], ax=axs[0])
# add annotations one by one with a loop
axs[0].title.set_text('best hits by iptm/plddt')
#axs[5].set_xlim([0.8, 1.0])
# add annotations one by one with a loop
for line in range(0,data_frame_top.shape[0]):
     axs[0].text(data_frame_top.iptm.iloc[line]+0.001, 
              data_frame_top.plddt.iloc[line],
              data_frame_top.run.iloc[line], 
              horizontalalignment='left',
              size='small', 
              color='black')

sns.scatterplot(data_frame_top['iptm'], data_frame_top['ptm'], ax=axs[1])
# add annotations one by one with a loop
axs[1].title.set_text('best hits by iptm/ptm')
#axs[5].set_xlim([0.8, 1.0])
# add annotations one by one with a loop
for line in range(0,data_frame_top.shape[0]):
     axs[1].text(data_frame_top.iptm.iloc[line]+0.001, 
              data_frame_top.ptm.iloc[line],
              data_frame_top.run.iloc[line], 
              horizontalalignment='left',
              size='small', 
              color='black')
     
plt.savefig('summary_best_scatter.png', dpi=300)

plt.clf()
sns.jointplot(data=data_frame, x="ptm", y="iptm")
plt.savefig('ptm_iptm_summary.png', dpi=300)


sorted_data_frame = data_frame.sort_values(by=['iptm'], ascending=False)
sorted_data_frame.to_csv('summary_sorted.csv', index=False)
