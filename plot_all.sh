#!/bin/bash
ls -d */ > list.txt
for i in `cat list.txt`
do 
	cp -v plot_AF_all.py $i/$i
	cp -v make_structure_figure.pml $i/$i
       	cd $i/$i 
	python plot_AF_all.py
	pymol -cq make_structure_figure.pml
       	pwd
       	cd -
done
