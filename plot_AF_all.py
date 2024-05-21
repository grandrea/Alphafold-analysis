import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
import pickle
from Bio import SeqIO
from itertools import accumulate
import json 
import glob
import numpy as np
import pandas as pd

# fasta file with input sequence must be in result directory
try:
    input_sequence_name = glob.glob('*fasta')[0]
except IndexError:
    #check if it is an alphafold 3 run
    try:
        input_sequence_name = glob.glob("*request*json")[0]
        alphafold3_json = True
    except IndexError:
        print("need .fasta file for alphafold2 run, .fasta or job request file for alphafold3")

#plot title as in fasta file
if alphafold3_json:
    print("parsing alphafold3 job request file\n")
    print("support currently limited to protein only")
    protein_names_for_title = input_sequence_name.replace("_job_request.json", "")
    print(protein_names_for_title)
    sequence_file = open(input_sequence_name)
    sequence_file_str = sequence_file.read()
    sequence_d = json.loads(sequence_file_str)
    sequence_list = []
    for idx, element in enumerate(sequence_d[0]["sequences"]):
        element_name = sequence_d[0]["name"] + "_chain_" + str(idx)
        element_sequence = sequence_d[0]["sequences"][idx]["proteinChain"]["sequence"]
        sequence_list.append((element_name, element_sequence))
    input_sequence = dict(sequence_list)
else:
    try:
        protein_names_for_title = input_sequence_name.replace('.fasta', '')
    except:
        protein_names_for_title = ''
    try:
        input_sequence = SeqIO.to_dict(SeqIO.parse(input_sequence_name, 'fasta'))
    # handle homomultimers with identical fasta headers
    except ValueError:
        with open(input_sequence_name) as original, open('corrected.fasta', 'w') as corrected:
            records = SeqIO.parse(original, 'fasta')
            for idx, record in enumerate(records):
                record.id = record.id + str(idx)
                print(record.id )
                SeqIO.write(record, corrected, 'fasta')
        corrected.close()
        input_sequence = SeqIO.to_dict(SeqIO.parse('corrected.fasta', 'fasta'))


# pickle files coming out of AF run for alphafold2
file_list = glob.glob('*_multimer*.pkl')
file_list = sorted(file_list)
alphafold3 = False
if len(file_list) < 1 : #assume no pickle files means an alphafold3 run
    if not alphafold3_json:
        print("guessing this is an AlphaFold3 directory\n")
        print("support currently limited to protein only")
    file_list = glob.glob('*full_data*.json')
    file_list = sorted(file_list)
    alphafold3 = True

# finished run will produce a .json file with the model ranking

if not alphafold3:
    try:
        json1_file = open('ranking_debug.json')
        json1_str = json1_file.read()
        model_stats = json.loads(json1_str)
        json1_file.close()
        model_ranking = model_stats.get('order')
    except FileNotFoundError:  # plot an incomplete run
        model_ranking = []

width_list = []
list_ranking = {}
statistics_list = []

if not alphafold3:
    # read pkl, match model number with rank using json file
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
else:
    for index, file_name in enumerate(file_list):
        stat_file = open(file_name)
        stat_file_str = stat_file.read()
        d1 = json.loads(stat_file_str)
        model_name = file_name.replace('.json', '')
        model_name = model_name.replace('_full_data_', '')

        summary_stat_name = file_name.replace("full_data", "summary_confidences")

        summary_file = open(summary_stat_name)
        summary_file_str = summary_file.read()
        d2 = json.loads(summary_file_str)

        d = dict(d1)
        d.update(d2)

        model_rank = index
        list_ranking[file_name] = [d, model_rank]
        width_list.append(4)
        statistics_list.append([model_name,
                                d.get('ptm'),
                                d.get('iptm'),
                                np.mean(d.get('atom_plddts')),
                                d.get('fraction_disordered'),
                                d.get("ranking_score")])

width_list.append(0.2)

# check protein lengths in input files for PAE plot
sequence_lengths = []

for item in input_sequence.keys():
    sequence = input_sequence.get(item)
    if not alphafold3_json:
        sequence_length = len(sequence.seq)
    else:
        sequence_length = len(sequence)
    sequence_lengths.append(sequence_length)
    

# figure out positions for black lines delimiting proteins in PAE plot
line_position_accumulated = accumulate(sequence_lengths)
line_positions = []
for item in line_position_accumulated:
    line_positions.append(item)


# write stats file--------
if not alphafold3:
    model_stats = pd.DataFrame(statistics_list,
                               columns = ['model', 'ptm', 'iptm', 'plddt', 'confrank'])
else:
    model_stats = pd.DataFrame(statistics_list,
                               columns=['model', 'ptm', 'iptm', 'plddt', 'fraction disordered', 'confrank'])

model_stats.to_csv('model_statistics.csv',  index=False)
# PAE plot -------------------

output_name = str('predicted_alignment_error.png')

palette = sns.diverging_palette(220, 20, as_cmap=True)

fig, axs = plt.subplots(ncols=len(file_list) + 1,
                        gridspec_kw=dict(width_ratios=width_list),
                        figsize = (4*len(file_list), 4))
fig.subplots_adjust(top=0.8)

for file_name in file_list:
    plot_number = list_ranking.get(file_name)[1]
    if not alphafold3:
        PAE = list_ranking.get(file_name)[0]['predicted_aligned_error']
        ipTM = list_ranking.get(file_name)[0]['iptm'].round(3)
    else:
        PAE = list_ranking.get(file_name)[0]['pae']
        ipTM = list_ranking.get(file_name)[0]['iptm']

    tick_range = [1] + list(range(500, len(PAE), 500))
    sns.heatmap(PAE,
                cmap=palette, 
                ax=axs[plot_number], 
                cbar=False, 
                vmin=0,
                vmax=30)
#    axs[plot_number].yaxis.set_major_locator(mticker.MaxNLocator(5))
#    axs[plot_number].xaxis.set_major_locator(mticker.MaxNLocator(5))
    if plot_number == 0:
        axs[plot_number].set_yticks(ticks = tick_range, labels = tick_range)
    else:
        axs[plot_number].set_yticks([], [])
    axs[plot_number].set_xticks(ticks=tick_range, labels=tick_range)
    axs[plot_number].title.set_text(str('model' + 
                                        str(list_ranking.get(file_name)[1]) + 
                                        '\n iptm: ' +
                                        str(ipTM)))
    
    # add black lines delimiting the two proteins
    for element in line_positions:
        axs[plot_number].vlines(element,
                                ymin=0,
                                ymax=len(PAE),
                                color='black',
                                linewidth=3)
        axs[plot_number].hlines(element,
                                xmin=0,
                                xmax=len(PAE),
                                color='black',
                                linewidth=3)
    
fig.colorbar(axs[0].collections[0], cax=axs[-1])
fig.suptitle(str('Predicted alignment error ' + protein_names_for_title))
plt.savefig(output_name)


# pLDDT plot------------------------------
if not alphafold3:
    output_name = str('pLDDT.png')

    # palette = sns.light_palette("#2ecc71", as_cmap=True, reverse=True)

    fig, axs = plt.subplots(ncols=len(file_list),
                            figsize = (6*len(file_list), 3))


    for file_name in file_list:
        plot_number = list_ranking.get(file_name)[1]
        PAE = list_ranking.get(file_name)[0]['plddt']
        tick_range = [1] + list(range(500, len(PAE), 500))
        ytick_range = list(range(0, 100, 10))
        try:
            axs[plot_number].plot(list(range(0,len(PAE), 1)), PAE, color='b')
            axs[plot_number].set_yticks(ticks=ytick_range, labels=ytick_range)
            axs[plot_number].set_xticks(ticks=tick_range, labels=tick_range)
            axs[plot_number].title.set_text(str('model' + str(list_ranking.get(file_name)[1])))
            for element in line_positions:
                axs[plot_number].vlines(element, ymin=0, ymax=100, color='black')
    #           #axs[plot_number].hlines(element, xmin = 0, xmax = len(PAE), color='black')
        except TypeError: #handle single model and single pkl file
            axs.plot(list(range(0,len(PAE), 1)), PAE, color='b')
            axs.set_yticks(ticks=ytick_range, labels=ytick_range)
            axs.set_xticks(ticks=tick_range, labels=tick_range)
            axs.title.set_text(str('model' + str(list_ranking.get(file_name)[1])))
            for element in line_positions:
                axs.vlines(element, ymin=0, ymax=100, color='black')

    fig.colorbar(axs[0].collections[0], cax=axs[-1])
    fig.suptitle(str('plddt ' + protein_names_for_title))
    plt.savefig(output_name)

    # MSA plot------------------------------

    def get_neff(msa, eff_cutoff=0.8): # eff_cutoff=0.62 for metapsicov
        if msa.ndim == 3: msa = msa.argmax(-1)
        # pairwise identity
        msa_sm = 1.0 - squareform(pdist(msa,"hamming"))
        # weight for each sequence
        msa_w = (msa_sm >= eff_cutoff).astype(float)
        msa_w = 1/np.sum(msa_w,-1)
        return msa_w.mean()




    if "features.pkl" in glob.glob("*"):
        features = pickle.load(open("features.pkl", "rb"))
        plt.close()
        plt.clf()
        ax = sns.heatmap(features["msa"])
        ax.set_xticks(ticks=tick_range, labels=tick_range)
        ax.set_xlabel("residue")
        ax.set_ylabel("sequences")
        plt.title("MSA information, logit units")
        for line_position in line_positions:
            plt.vlines(line_position,
                       ymin=0,
                       ymax=len(features["msa"][:,1]),
                       color="black",
                       linewidth=3)
        plt.xticks(tick_range)
        plt.subplots_adjust(bottom=0.15)
        plt.savefig("MSA.png")


