# Alphafold-analysis


Generate statistics plots for AlphaFold2 and AlphaFold3 protein multimer runs.


Requires numpy, seaborn, pandas,  matplotlib, biopython. It may require jax to unpack .pkl files produced by some versions of AlphaFold2. Pymol is optional and used by make_structure_figure.pml.

Collection of scripts to survey a directory of alphafold2 multimer and alphafold3 runs. For alphafold3, download the run from alphafold server. Only protein-protein runs are supported for alphafold3. Currently, the wrapper scripts only work with alphafold2.

Put the .fasta file with the sequence submitted to alphafold in the results folder. For alphafold3, sequence is parsed directly from the job request json file coming from alphafold server.

To analyse a single directory, no wrapper scripts:

```
python plot_AF_all.py
```

will make PAE, pLDDT and MSA coverage plots.

![PAE plot](https://i.imgur.com/f41BenC.png)

```
pymol -cq make_structure_figure.pml
```
will make a figure of the structure colored by pLDDT, one colored by chains and one showing hydrogen bonds. Edit the first line to load a .cif file from alphafold3. 
This is at the moment coded for dimers.


#### Split AlphaFold 2.x and 3.x models by confidence. 

This is a good remedy for AF3 "hallucinating" helices from disordered regions- it essentially removes from the model all low confidence regions. Works on both pdb and cif files.

Run with

    python alphacrop.py --filename mystructure.pdb

or mystructure.cif. If needed, explicitly specify file format with the --file_format flag

Confidence defined by Deepmind, each input into multiple files made up of:
- only very high confidence regions (plDDT > 90)
- only confident regions and above (plDDT > 70)
- only low concfidence regions and above, discarding very low confidence (plDDT > 50)

Essentially splits the file by b-factor column, generating 3 files.

Also available as colab here  https://colab.research.google.com/drive/1kfF_Kkozi5ox8mVebqKgvc7UEWdXWpnf?usp=sharing

#### For Alphafold2 only

If instead you have a whole directory of runs, make sure each run has a name like RAB-RAK (for dimers), or similar, with the protein names, and identical subfolder name.
The script  will go in each directory and make pLDDT and PAE plots delimited by sequence, provided the results directory contains a .fasta of the sequence.

first make all summary files and plot all PAE, pLDDT by

```
./plot_all.sh
```
then make a summary file

```
./make_summary.sh
```
then create a summary plot

```
python alphafold_summary.py
```

