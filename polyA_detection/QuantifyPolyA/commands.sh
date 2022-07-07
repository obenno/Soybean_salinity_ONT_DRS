#! /usr/bin/env bash

## Generate QuantifyPolyA input bed file
for i in ../{Leaf,Root}*R{1,2,3}.PAS.out
do
    filename=$(basename $i)
    awk 'NR>1{print $3"\t"$4"\t"$5}' $i | sort | uniq -c |
        awk '{print $2"\t"$4"\t"$3"\t"$1}' |
        sort -k 1,1 -k2,2 -k 3,3n > ${filename%%.out}".bed"
done

## Identify polyA sites by QuantifyPolyA
Rscript QuantifyPolyA.R

## Aggregate polyA information for all the PACs
for i in ../*R{1,2,3}.PAS.out
do
    sampleName=$(basename $i)
    sampleName=${sampleName%%.PAS.out}
    awk 'NR>1{print $3"\t"$4-1"\t"$4"\t"$1"\t"$2"\t"$5}' $i |
        bedtools intersect -a QuantifyPolyA_polyA_sites.bed -b - -wo |
        cut -f 1-6,11| sort -k 4,4 -k 7,7n |
        bedtools groupby -g 1,2,3,4,5,6 -c 7 -o mean > $sampleName".polyALen.tsv"
done

## Plot polyA length boxplot
## Categories proximal and distal PACs
Rscript plot_polyAlen_boxplot.R
