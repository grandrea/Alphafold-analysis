# Alphafold-analysis

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

