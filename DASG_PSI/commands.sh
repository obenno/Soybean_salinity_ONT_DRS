#! /usr/bin/env bash

## Aggregate AS and aPAS events expressions
cat ../alternative_splicing/AS_events_tpm.tsv ../polyA_detection/aPAS_events_cpm.tsv > AS_aPAS_event_expr.tsv

## Stats for delta PSI value
## This script calculate PSI, deltaPSI and conduct a fisher exact test
## dAS events (including dPAS) were selected by:
## |deltaPSI| >= 0.1 and FDR <= 0.05
## outputs: DASG.lst, dAS_eventCount.tsv, dAS_geneCount.tsv, PSI.pdf
Rscript calc_dAS_Pvalue.R

## Gene ref gene list
awk '$3=="transcript"{print $10"\t"$12}' ../../stringtie_pipeline/stringtie.annotated.gtf | sed 's/\"//g; s/;//g' | awk 'NR==FNR{a[$1]=$2}NR>FNR{if($1 in a){print $2"\t"$1"\t"a[$1]}}' ../../stringtie_pipeline/correspondingList - | awk '{split($3, tmp, "."); print $1"\t"tmp[1]"."tmp[2]}'| sort -u > known_genes.lst

## copy DEG.lst from illumina DE result
cp ../../expressions/illumina/DEG.lst ./

## Extract DEG and DASG from refGene list
awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print $2}}' DEG.lst known_genes.lst > DEG_refGene.lst
awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print $2}}' DASG.lst known_genes.lst > DASG_refGene.lst

## Plot DEG DASG venn diagram
Rscript vennPlot.R -i DEG.lst,DASG.lst -n DEG,DASG -o DEG_vs_DASG.venn.pdf

## Count differential events with at least one isoform completely from novel transcripts
awk '$3=="transcript"{print $12"\t"$1"\t"$7}' ../stringtie_pipeline/stringtie.annotated.gtf | sed 's/\"//g; s/;//' | sort -u | awk 'NR==FNR{a[$1]=$2;b[$1]=$3}NR>FNR{print $1"\t"$2"\t"a[$1]":"$3"-"$4":"b[$1]"\t"a[$1]":"$5"-"$6":"b[$1]}' -  ../alternative_splicing/AS_events_withNovelIsoform.tsv | awk 'NR==FNR{a[$0]}NR>FNR{if($1"\t"$2"\t"$3"\t"$4 in a){print}}' - dAS_data.tsv | awk '{print $2}'| sort | uniq -c
