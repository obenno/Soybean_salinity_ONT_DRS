#! /usr/bin/env bash

## Transdecoder commands
~/Tools/TransDecoder/util/gtf_to_alignment_gff3.pl stringtie.annotated.gtf > stringtie.transcripts.gff3

~/Tools/TransDecoder/TransDecoder.LongOrfs -t stringtie.annotated.transcript.fa

~/Tools/TransDecoder/TransDecoder.Predict -t stringtie.annotated.transcript.fa

~/Tools/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl stringtie.annotated.transcript.fa.transdecoder.gff3 stringtie.transcripts.gff3 stringtie.annotated.transcript.fa > stringtie.annotated.transcript.fa.transdecoder.genome.gff3

## Select the longest CDS and add to bed20 format
awk '$3=="CDS"{split($9, attr, ";"); id=substr(attr[2],8); split(id, tmp, "."); print tmp[1]"\t"tmp[2]"\t"$5-$4+1}' stringtie.annotated.transcript.fa.transdecoder.gff3 | sort -k 1,1 -k 3,3rn | bedtools groupby -g 1 -c 2 -o collapse | awk '{split($2, tmp, ","); print $1"."tmp[1]}'| awk 'NR==FNR{a[$1]}NR>FNR && $3=="CDS"{split($9, attr, ";"); id=substr(attr[2],8); if(id in a){print id"\t"$4-1"\t"$5}}' - stringtie.annotated.transcript.fa.transdecoder.genome.gff3 | sort -k1,1 -k 2,2n -k 3,3n |bedtools groupby -g 1 -c 2,3 -o min,max | awk '{split($1, tmp, "."); print tmp[1]"\t"$2"\t"$3}' |awk 'BEGIN{FS="\t"; OFS="\t"}NR==FNR{s[$1]=$2; e[$1]=$3}NR>FNR{if($4 in s){$7=s[$4]; $8=e[$4]}; print}' - ../stringtie_pipeline/stringtie.count_5.colored.bed20 > stringtie.count_5.colored.withCDS.bed20

## Correct some delimiter issue
sed -i 's/\t\t/\t/' stringtie.count_5.colored.withCDS.bed20

## Select the longest CDS and add to gtf format
## Gtf is used to make txdb object in subsequent analysis
awk '$3=="CDS"{split($9, attr, ";"); id=substr(attr[2],8); split(id, tmp, "."); print tmp[1]"\t"tmp[2]"\t"$5-$4+1}' stringtie.annotated.transcript.fa.transdecoder.gff3 | sort -k 1,1 -k 3,3rn | bedtools groupby -g 1 -c 2 -o collapse | awk '{split($2, tmp, ","); print $1"."tmp[1]}'| awk 'NR==FNR{a[$1]}NR>FNR && $3=="CDS"{split($9, attr, ";"); id=substr(attr[2],8); if(id in a){print id"\t"$1"\t"$4"\t"$5"\t"$7"\t"$8}}' - stringtie.annotated.transcript.fa.transdecoder.genome.gff3 | awk '{split($1, tmp, "."); print tmp[1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print}}' ../stat/stringtie.transcript.counts_more_than_five.lst -|awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]"\t"$0}' ../polyA_detection/trans_gene.tsv - | awk '{print $3"\ttransdecoder\tCDS\t"$4"\t"$5"\t.\t"$6"\t"$7"\ttranscript_id \""$2"\"; gene_id \""$1"\";"}' | cat - ../stringtie_pipeline/stringtie.count_5.gtf - | sort -k 1,1 -k 4,4n > stringtie.count_5.withCDS.gtf

## Replace unknown strand infor with "+"
awk '$7=="."{$7="+"}{print}' FS="\t" OFS="\t" stringtie.count_5.withCDS.gtf > stringtie.count_5.withCDS.strandCorrected.gtf
