# Alphafold-analysis

Requires numpy, seaborn, pandas, numpy, matplotlb, biopython

Collection of scripts to survey a directory of alphafold multimer runs.

it will go in each directory and make pLDDT and PAE plots delimited by sequence, provided the results directory contains a .fasta of the sequence.

first make all summary files and plot all PAE, pLDDT by

./plot_all.sh

then make a summary file

./make_summary.sh

then create a summary plot

python alphafold_summary.py


