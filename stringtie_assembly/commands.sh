#! /usr/bin/env bash

## StringTie commands
stringtie -G ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -p 40 -L -o singleLib/C08leaf_0h.stringtie.gtf ../minimap2_mapping/C08leaf_0h.minimap2.bam
stringtie -G ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -p 40 -L -o singleLib/C08leaf_1h.stringtie.gtf ../minimap2_mapping/C08leaf_1h.minimap2.bam
stringtie -G ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -p 40 -L -o singleLib/C08root_1h.stringtie.gtf ../minimap2_mapping/C08root_1h.minimap2.bam
stringtie -G ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -p 40 -L -o singleLib/C08root_0h.stringtie.gtf ../minimap2_mapping/C08root_0h.minimap2.bam
stringtie --merge -G ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -F 0 -T 0 -i -o merged.gtf *.gtf
## Collect isoforms not included in stringtie merge
gffcompare -r merged.gtf C08*.gtf
## Comapre and annotate assembled transcripts with annotation
gffcompare -r ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 gffcmp.combined.gtf -o stringtie
## Generate transcripts fasta file
gffread -g ~/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa -w stringtie.annotated.transcript.fa stringtie.annotated.gtf

## Filter stringtie transcript by reads count >= 5
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 ../stat/nanoReads_vs_StringTieTranscriptome.minimap2.bam | awk '{print $1"\t"$3}' | sort -k 2,2 | bedtools groupby -g 2 -c 1 -o count | awk '$2>=5{print $1"\t"$2}' > stringtie.transcript.counts_more_than_five.lst

awk 'NR==FNR{a["\""$1"\";"]}NR>FNR{if($10 in a){print}}' ../stat/stringtie.transcript.counts_more_than_five.lst stringtie.annotated.gtf > stringtie.count_5.gtf

## Prepare colored bed file for jbrowse and IGV

## Generate correspondingList for annotated transcripts
## For non filtered transcripts gtf
awk -F"\t" '$3=="transcript"{split($9,tmp,"; "); for(i=1;i<=length(tmp);i++){split(tmp[i],k," "); a[k[1]]=k[2]};if(a["class_code"]=="\"=\""){print a["transcript_id"]"\t"a["cmp_ref"]}}' stringtie.annotated.gtf |  sed 's/\"//g' > correspondingList

## For filtered transcripts gtf
gtfToGenePred -genePredExt stringtie.count_5.gtf stringtie.count_5.gp
genePredToBigGenePred stringtie.count_5.gp stdout | sort -k1,1 -k2,2n > stringtie.count_5.bed20
awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1]}NR>FNR{if(!($4 in a)){if($6=="+"){$9="239,138,98"}else{$9="103,169,207"}}; print}' correspondingList stringtie.count_5.bed20 > stringtie.count_5.colored.bed20

