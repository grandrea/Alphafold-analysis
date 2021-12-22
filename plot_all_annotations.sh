#!/bin/bash
ls -d */ > list.txt
for i in `cat list.txt`
do 
	cp -v plot_AF_all.py $i/$i
       	cd $i/$i 
	FILE=./predicted_alignment_error.png
	if [ ! -f "$FILE" ]; then
		echo "plotting in `pwd`"
		python plot_AF_all.py
		signalp -fasta *fasta -format long
		cat *.fasta | /home/andrea/software/tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64 -f /home/andrea/software/tmhmm-2.0c/lib/TMHMM2.0.options -modelfile /home/andrea/software/tmhmm-2.0c/lib/TMHMM2.0.model > transmbembrane.txt
	fi
       	cd -
done
