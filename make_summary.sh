#!/bin/bash

grep model ./*/*/*.csv > summary.txt

grep -v "model,ptm,iptm,plddt,confrank" summary.txt > summary.csv

rm summary.txt

sed -i 's/model_statistics.csv:/,/g' summary.csv
sed -i 's/\/,/,/g' summary.csv
sed -i 's/^.*\///g' summary.csv

echo "run,model,ptm,iptm,plddt,confrank" > a
cat a summary.csv > b
rm a
mv b summary.csv


