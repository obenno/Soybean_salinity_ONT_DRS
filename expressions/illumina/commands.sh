#! /usr/bin/env bash

## Raw reads QC steps
ls *_1.fastq.gz | awk '{sample=substr($1,1,index($1,"_1.fastq.gz")-1); print "trimmomatic PE -threads 30 -phred33 "$1" "sample"_2.fastq.gz "sample"_1.trimmed.fastq.gz "sample"_1.unpaired.fastq.gz "sample"_2.trimmed.fastq.gz "sample"_2.unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 MINLEN:36"}' > QC.sh
## perform QC
bash QC.sh

## Get decoys.txt
fasta_formatter -t -i ~/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa | awk '{print $1}' > decoys.txt
## Concatenate transcripts fasta and genome fasta
cat stringtie.annotated.transcript.fa ~/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa | gzip > gentrome.fa.gz
## Build salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -p 20 -i salmon_index --keepDuplicates
## Generate salmon script
ls *_1.trimmed.fastq.gz | awk '{sample=substr($1, 1, index($1, "_1.trimmed.fastq.gz")-1); print "salmon quant -i salmon_index -l A -1 "$1" -2 "sample"_2.trimmed.fastq.gz -p 40 --validateMappings -o "sample}' > salmon_commands.sh
## Run salmon
bash salmon_commands.sh

## Prepare DE counts
## For DEG/DEI analysis following salmon, info files need to be prepared:
## transcript to gene infor: tx2gene.tsv
awk '$3=="transcript"{print $10"\t"$12}' ../../stringtie_pipeline/stringtie.annotated.gtf | sed 's/\"//g; s/;//g' | awk 'BEGIN{print "TXNAME\tGENEID"}NR==FNR{a[$1]=$2}NR>FNR && FNR>1{print $1"\t"a[$1]}' - SRR6669437/quant.sf | awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR==1{print $0"\tRef_transcript"}NR>FNR && FNR>1{if($1 in a){print $0"\t"a[$1]}else{print $0"\tNA"}}' ../../stringtie_pipeline/correspondingList - > tx2gene.tsv

## Perform DEG and DEI analysis
Rscript calulateDE.R

## Calculate DEG and DEIG numbers
awk 'FNR>1{print $1}' C08*gene*.out | sort -u > DEG.lst
awk 'FNR>1{print $1}' C08*transcript*.out | sort -u | awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' tx2gene.tsv - | sort -u > DEIG.lst
## Plot venn diagram
Rscript vennPlot.R -i DEG.lst,DEIG.lst -n DEGs,DEIGs -o DEG_vs_DEIG.venn.pdf
