import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
import pickle
from Bio import SeqIO
from itertools import accumulate
import json 
import glob
import numpy as np
import pandas as pd

#fasta file with input sequence must be in result directory
input_sequence_name  = glob.glob('*fasta')[0]

#plot title as in fasta file
try:
    protein_names_for_titile = input_sequence_name.replace('.fasta', '')
except:
    protein_names_for_titile = ''


try:
    input_sequence = SeqIO.to_dict(SeqIO.parse(input_sequence_name, 'fasta'))
#handle homomultimers with identical fasta headers
except ValueError:
    with open(input_sequence_name) as original, open('corrected.fasta', 'w') as corrected:
        records = SeqIO.parse(original, 'fasta')
        for idx, record in enumerate(records):
            record.id = record.id + str(idx)
            print(record.id )
            SeqIO.write(record, corrected, 'fasta')
    corrected.close()
    input_sequence = SeqIO.to_dict(SeqIO.parse('corrected.fasta', 'fasta'))

#pickle files coming out of AF run
file_list = glob.glob('*_multimer.pkl')
file_list = sorted(file_list)

#finished run will produce a .json file with the model ranking 

try:
    json1_file = open('ranking_debug.json')
    json1_str = json1_file.read()
    model_stats = json.loads(json1_str)
    json1_file.close()
    model_ranking = model_stats.get('order')
except FileNotFoundError: #plot an incomplete run 
    model_ranking = []



PAE_list = []

width_list = []

list_ranking = {}

statistics_list = []

#read pkl, match model number with rank using json file
for index, file_name in enumerate(file_list):
    d = pickle.load(open(file_name,'rb'))
    model_name = file_name.replace('.pkl', '')
    model_name = model_name.replace('result_', '')
    try:
        model_rank = model_ranking.index(model_name)
    except ValueError:
        model_rank = index
    list_ranking[file_name] = [d, model_rank]
    width_list.append(4)
    statistics_list.append([model_name, 
                            d.get('ptm'), 
                            d.get('iptm'),
                            np.mean(d.get('plddt')),
                            d.get('ranking_confidence')])
    PAE_list.append(d)

width_list.append(0.2)

#check protein lengths in input files for heatmap
sequence_lengths = []

for item in input_sequence.keys():
    sequence = input_sequence.get(item)
    sequence_length = len(sequence.seq)
    sequence_lengths.append(sequence_length)
    


line_position_accumulated = accumulate(sequence_lengths)
line_positions = []
for item in line_position_accumulated:
    line_positions.append(item)


#write stats file--------
model_stats = pd.DataFrame(statistics_list, 
             columns = ['model', 'ptm', 'iptm', 'plddt', 'confrank'])
model_stats.to_csv('model_statistics.csv',  index=False )
#PAE plot -------------------

output_name = str('predicted_alignment_error.png')

palette = sns.diverging_palette(220, 20, as_cmap=True)

fig, axs = plt.subplots(ncols=len(file_list) + 1,
                        gridspec_kw=dict(width_ratios=width_list),
                        figsize = (4*len(file_list), 4))


for file_name in file_list:
    plot_number = list_ranking.get(file_name)[1]
    PAE = list_ranking.get(file_name)[0]['predicted_aligned_error']
    tick_range = [1] + list(range(500, len(PAE), 500))
    sns.heatmap(PAE, cmap=palette, ax=axs[plot_number], cbar=False, vmin= 0, vmax = 30)
#    axs[plot_number].yaxis.set_major_locator(mticker.MaxNLocator(5))
#    axs[plot_number].xaxis.set_major_locator(mticker.MaxNLocator(5))
    if plot_number == 0:
        axs[plot_number].set_yticks(ticks = tick_range, labels = tick_range)
    else:
        axs[plot_number].set_yticks([], [])
    axs[plot_number].set_xticks(ticks = tick_range, labels = tick_range)
    axs[plot_number].title.set_text(str('model' + str(list_ranking.get(file_name)[1])))
    for element in line_positions:
        axs[plot_number].vlines(element, ymin = 0, ymax = len(PAE), color='black')
        axs[plot_number].hlines(element, xmin = 0, xmax = len(PAE), color='black')
    
fig.colorbar(axs[0].collections[0], cax=axs[-1])
fig.suptitle(str('Predicted alignment error ' + protein_names_for_titile)) 
plt.savefig(output_name)


#pLDDT plot------------------------------

output_name = str('pLDDT.png')


#palette = sns.light_palette("#2ecc71", as_cmap=True, reverse=True)

fig, axs = plt.subplots(ncols=len(file_list),
                        figsize = (6*len(file_list), 3))


for file_name in file_list:
    plot_number = list_ranking.get(file_name)[1]
    PAE = list_ranking.get(file_name)[0]['plddt']
    tick_range = [1] + list(range(500, len(PAE), 500))
    axs[plot_number].plot(list(range(0,len(PAE), 1)), PAE, color='b')
    if plot_number == 0:
        axs[plot_number].set_yticks(ticks = tick_range, labels = tick_range)
    else:
       axs[plot_number].set_yticks([], [])
    axs[plot_number].set_xticks(ticks = tick_range, labels = tick_range)
    axs[plot_number].title.set_text(str('model' + str(list_ranking.get(file_name)[1])))
    for element in line_positions:
        axs[plot_number].vlines(element, ymin = 0, ymax = 100, color='black')
#         #axs[plot_number].hlines(element, xmin = 0, xmax = len(PAE), color='black')
    
#fig.colorbar(axs[0].collections[0], cax=axs[-1])
fig.suptitle(str('plddt ' + protein_names_for_titile)) 
plt.savefig(output_name)

