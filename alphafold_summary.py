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

#single model with best iptm
data_frame_onemodel_iptm = data_frame.sort_values(by=['iptm'], ascending=False)
data_frame_onemodel_iptm = data_frame_onemodel_iptm.drop_duplicates(subset = ["run"])

plt.clf()
fig, axs = plt.subplots(ncols=3,
                        figsize = (16, 3))

sns.histplot(data=data_frame, x='ptm', bins=20, ax=axs[0])
axs[0].title.set_text('pTM score')


sns.histplot(data=data_frame, x='iptm', bins=20, ax=axs[1])
axs[1].title.set_text('ipTM score')

sns.histplot(data=data_frame,x='plddt', bins=20, ax=axs[2])
axs[2].title.set_text('plddt')


plt.suptitle('Alphfold run summary best iptm structures only')
plt.savefig('summary_histograms.png', dpi=300)

plt.clf()
fig, axs = plt.subplots(ncols=3,
                        figsize = (16, 3))
sns.scatterplot(data=data_frame_onemodel_iptm, x='ptm', y='plddt', ax=axs[0])
axs[0].title.set_text('plddt vs pTM')
axs[0].set_xlim(0,1)
axs[0].set_ylim(0,100)

sns.scatterplot(data=data_frame_onemodel_iptm, x='iptm', y='plddt', ax=axs[1])
axs[1].title.set_text('plddt vs ipTM')
axs[1].set_xlim(0,1)
axs[1].set_ylim(0,100)

sns.scatterplot(data=data_frame_onemodel_iptm, x='ptm', y='iptm', ax=axs[2])
axs[2].title.set_text('ipTM vs pTM')
axs[2].set_xlim(0,1)
axs[2].set_ylim(0,1)


plt.suptitle('Alphfold run summary')
plt.savefig('summary_scatter.png', dpi=300)

plt.clf()
fig, axs = plt.subplots(ncols=2,
                        figsize = (12, 2))

sns.scatterplot(data=data_frame_top, x='iptm', y='plddt', ax=axs[0])
# add annotations one by one with a loop
axs[0].title.set_text('best hits by iptm/plddt')
for line in range(0,data_frame_top.shape[0]):
     axs[0].text(data_frame_top.iptm.iloc[line]+0.001, 
              data_frame_top.plddt.iloc[line],
              data_frame_top.run.iloc[line], 
              horizontalalignment='left',
              size='small', 
              color='black')

sns.scatterplot(data=data_frame_top, x='ptm', y='iptm', ax=axs[1])
# add annotations one by one with a loop
axs[1].title.set_text('best hits by iptm/ptm')
for line in range(0,data_frame_top.shape[0]):
     axs[1].text(data_frame_top.ptm.iloc[line]+0.001, 
              data_frame_top.iptm.iloc[line],
              data_frame_top.run.iloc[line], 
              horizontalalignment='left',
              size='small', 
              color='black')
     
plt.savefig('summary_best_scatter.png', dpi=300)

plt.clf()
sns.jointplot(data=data_frame_onemodel_iptm, x="ptm", y="iptm", xlim = (0,1),ylim = (0,1))
plt.savefig('ptm_iptm_summary.png', dpi=300)


sorted_data_frame = data_frame.sort_values(by=['iptm'], ascending=False)
sorted_data_frame.to_csv('summary_sorted.csv', index=False)
