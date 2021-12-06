#!/bin/bash
ls -d */ > list.txt
for i in `cat list.txt`
do 
	cp -v plot_AF_all.py $i/$i
       	cd $i/$i 
	python plot_AF_all.py
       	pwd
       	cd -
done
