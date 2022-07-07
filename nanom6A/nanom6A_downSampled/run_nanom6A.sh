#! /usr/bin/env bash

for i in `ls | awk '$1~/Leaf|Root/'`
do
    predict_sites --cpu 40 \
                  -i $i/result \
                  -o $i/result_final \
                  -r /home/ubuntu/salinity_suppl_analysis/analysis/stringtie/stringtie.annotated.gene.bed \
                  -g ~/data/database/Gmax_508_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa \
                  --proba 0.5 \
                  --support 10
done
