#! /usr/bin/env bash

## Commands for m6A sites and ratio stats

## Generate whole genome m6A sites (without depth filtration)
function tidy_m6A_coordinates {
    sort --parallel=8 -S 10% -k 3,3 -k 1,1 -k 2,2n -k 4,4 ../$1/result_final/genome_abandance.0.5.bed |
        cut -f 1,2,6 |
        sort -u --parallel=8 -S 10% | sort -k 1,1 -k 2,2n --parallel=8 -S 10% > m6A_wholeGenomeSites.$1.tsv
}

for i in `ls ../ | awk '$1~/Leaf|Root/'`;
do
    sampleName=$i
    tidy_m6A_coordinates $sampleName
done

cat m6A_wholeGenomeSites.*.tsv | sort -u --parallel=8 -S 10% | sort -k 1,1 -k 2,2n --parallel=8 -S 10% > m6A_wholeGenomeSites.tsv

## Plot m6A sites ratio distribution, m6A sites were also extracted (depth >=20)
function tidy_m6A_info {
    awk '{split($1, tmp, "|"); gene=tmp[1]; chr=tmp[2]; for(i=2;i<=NF;i++){split($i, d, "|"); pos=d[1]; methyl_count=d[2]; total_count=d[3]; ratio=d[4]; print gene"\t"chr"\t"pos"\t"ratio"\t"methyl_count"\t"total_count}}' ../$1/result_final/ratio.0.5.tsv |
        awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2"\t"$3]"\t"$4"\t"$5"\t"$6}' m6A_wholeGenomeSites.tsv - > m6A_ratio.$1.tsv
}

for i in `ls ../ | awk '$1~/Leaf|Root/'`;
do
    sampleName=$i
    tidy_m6A_info $sampleName
done


#### aggregate m6A sites
Rscript find_rep_common.R Leaf_0h_R1 Leaf_0h_R2 Leaf_0h_R3 Leaf_0h.ratioCommon.tsv
Rscript find_rep_common.R Leaf_1h_R1 Leaf_1h_R2 Leaf_1h_R3 Leaf_1h.ratioCommon.tsv
Rscript find_rep_common.R Root_0h_R1 Root_0h_R2 Root_0h_R3 Root_0h.ratioCommon.tsv
Rscript find_rep_common.R Root_1h_R1 Root_1h_R2 Root_1h_R3 Root_1h.ratioCommon.tsv

## Plot bar plots for calculating m6A sites
## and generate m6A site lists
Rscript plot_m6A_sites_bar.R

## Plot venn plot of nodule and sr
Rscript vennPlot.R -i Leaf_0h_sites.lst,Leaf_1h_sites.lst,Root_0h_sites.lst,Root_1h_sites.lst \
        -n Leaf_0h,Leaf_1h,Root_0h,Root_1h -o m6A_sites.venn.pdf \
        --width 2.5 --height 2.5
rm VennDiagram*.log

## Plot replicate correlation, boxplot and scatter plot
Rscript plot_repCor.R

## ## Generate gene TPM table
awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR>1{print a[$1]"\t"$0}' ../../polyA_detection/trans_gene.tsv ../../expressions/Nano_TPM.long.tsv > Nanopore_salmonminimap2_TPM.long.tsv
## ## Plot ratio vs gene expression
Rscript plot_m6A_ratio_vs_geneExpr.R

## Generate m6A sites bed files
## and plot m6A distribution in gene regions (5'UTR, CDS, 3'UTR)
awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' Leaf_0h_sites.lst > Leaf_0h_sites.bed
awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' Leaf_1h_sites.lst > Leaf_1h_sites.bed
awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' Root_0h_sites.lst > Root_0h_sites.bed
awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' Root_1h_sites.lst > Root_1h_sites.bed

## Plot m6A position density across gene of the four samples together in one plot
Rscript plot_m6A_geneFeatures_density.R -i ../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf -s Leaf_0h_sites.bed,Leaf_1h_sites.bed,Root_0h_sites.bed,Root_1h_sites.bed -n Leaf_0h,Leaf_1h,Root_0h,Root_1h -t 30 -o m6A.GenomicFeatures_density.allSamples.pdf

## Plot m6A position density of Leaf_0h_specific, Leaf_common and Leaf_1h_specific
awk 'NR==FNR{a[$1]}NR>FNR{if(!($1 in a)){print}}' Leaf_1h_sites.lst Leaf_0h_sites.lst |
    awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Leaf_1hvs0h.0h_specific.bed

awk 'NR==FNR{a[$1]}NR>FNR{if(!($1 in a)){print}}' Leaf_0h_sites.lst Leaf_1h_sites.lst |
    awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Leaf_1hvs0h.1h_specific.bed

awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print}}' Leaf_1h_sites.lst Leaf_0h_sites.lst | awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Leaf_1hvs0h.common.bed

Rscript plot_m6A_geneFeatures_density.R -i ../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf -s Leaf_1hvs0h.0h_specific.bed,Leaf_1hvs0h.common.bed,Leaf_1hvs0h.1h_specific.bed -n Leaf_0h_specific,Leaf_common,Leaf_1h_specific -t 30 -o m6A.GenomicFeatures_density.leafSpecific.pdf

## Plot m6A position density of Root_0h_specific, Root_common and Root_1h_specific
awk 'NR==FNR{a[$1]}NR>FNR{if(!($1 in a)){print}}' Root_1h_sites.lst Root_0h_sites.lst |
    awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Root_1hvs0h.0h_specific.bed

awk 'NR==FNR{a[$1]}NR>FNR{if(!($1 in a)){print}}' Root_0h_sites.lst Root_1h_sites.lst |
    awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Root_1hvs0h.1h_specific.bed

awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print}}' Root_1h_sites.lst Root_0h_sites.lst | awk '{split($1, tmp, "::");  print tmp[2]"\t"tmp[3]-1"\t"tmp[3]"\t"$1}' > Root_1hvs0h.common.bed

Rscript plot_m6A_geneFeatures_density.R -i ../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf -s Root_1hvs0h.0h_specific.bed,Root_1hvs0h.common.bed,Root_1hvs0h.1h_specific.bed -n Root_0h_specific,Root_common,Root_1h_specific -t 30 -o m6A.GenomicFeatures_density.rootSpecific.pdf

