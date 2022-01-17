# Alphafold-analysis

Requires numpy, seaborn, pandas,  matplotlb, biopython, pymol

Collection of scripts to survey a directory of alphafold multimer runs.

To analyse a single multimer run, no wrapper scripts:

```
python plot_AF_all.py
```

will make the PAE and pLDDT plots.

![PAE plot](https://imgur.com/a/33xH3Bx)

```
pymol -cq make_structure_figure.pml
```
will make a figure of the structure colored by pLDDT, one colored by chains and one showing hydrogen bonds.
This is at the moment coded for dimers.

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

