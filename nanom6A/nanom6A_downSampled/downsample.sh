#! /usr/bin/env bash

for i in `ls | awk '$1~/Leaf|Root/'`;
do
    shuf -n 976779 /home/ubuntu/salinity_suppl_analysis/analysis/nanom6A/$i/files.txt > $i/files.txt
done
