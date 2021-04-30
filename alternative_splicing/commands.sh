#! /usr/bin/env bash

## These commands were used to get PSI result from AS events
## other than aPAS

## Filter stringtie result by transcript count >= 5
awk 'NR==FNR{a["\""$1"\";"]}NR>FNR{if($10 in a){print}}' ../stat/stringtie.transcript.counts_more_than_five.lst stringtie.annotated.gtf > stringtie.count_5.gtf

## Detect and extract AS events from stringtie gtf
python extract_AS_events.py -i stringtie.count_5.gtf -o stringtie.count_5.AS.out

## Ensure at least one inclusion isoform and one exclusion isoform has max TPM >= 1
awk 'NR==FNR{a[$1]}NR>FNR{split($7, inclusion, ","); split($8, exclusion, ","); m=0; n=0; for(i=1;i<=length(inclusion);i++){if(inclusion[i] in a){m=1}}; for(i=1;i<=length(exclusion);i++){if(exclusion[i] in a){n=1}}; if(m==1 && n==1){print}}' ../expressions/Nanopore/salmonminimap2/Nano_maxTPM1.lst stringtie.count_5.AS.out > stringtie.count_5.maxTPM_1.AS.out

## Format AS events inclusion and exclusion expression table
## for subsequent fisher test
awk '$3=="transcript"{print $1"\t"$7"\t"$12}' stringtie.count_5.gtf | sed 's/\"//g; s/;//'|sort -u | awk 'NR==FNR{a[$3]=$1"\t"$2}NR>FNR{split(a[$1], tmp,"\t"); chr=tmp[1]; strand=tmp[2]; print $1"\t"$2"\t"chr":"$3"-"$4":"strand"\t"chr":"$5"-"$6":"strand"\t"$7"\t"$8}' - <( awk -f format_AS_events_expr.awk ../expressions/Nanopore/salmonminimap2/Nano_TPM.long.tsv stringtie.count_5.maxTPM_1.AS.out) > AS_events_tpm.tsv

## Extract splice sites from annotation, stringtie result
Rscript extract_intronRegion.R -i ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -g ~/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa -o Gmax_a2v1_annotation_introns.tsv

Rscript extract_intronRegion.R -i stringtie.count_5.strandEnforced.gtf -o stringtie.count_5.introns.tsv -g ~/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa

## Perform splice sites statistics, grouped by "Canonical" and "Non-Canonical"
## Check how many were annotated and how many were confirmed by short reads
awk 'FILENAME==ARGV[1]{anno[$1"\t"$2"\t"$3]}FILENAME==ARGV[2]{short[$1"\t"$2"\t"$3]}FILENAME==ARGV[3]{if($1"\t"$2"\t"$3 in anno){print $0"\tAnnotated"}else if($1"\t"$2"\t"$3 in short){print $0"\tShortReads"}else{print $0"\tNovel"}}' Gmax_a2v1_annotation_introns.tsv ../expressions/shortReads_STAR_mapping/aggregated_STAR_spliceSite.tsv <(awk '{print $1"\t"$2"\t"$3"\t"$NF}' stringtie.count_5.introns.tsv) > splice_sites_stat.input.tsv

## Plot bar and output count table
## output spliceSite_stat.out and spliceSite_annotated.bar.pdf
Rscript plot_spliceSite_bar.R

## Count how many events containing one isofrom completely from novel transcripts
awk 'NR==FNR{k[$1]=$2}NR>FNR&&FNR>1{split($7, inclusion, ","); split($8, exclusion, ","); a=0; b=0; for(i=1;i<=length(inclusion);i++){if(inclusion[i] in k){a=1}}; for(i=1;i<=length(exclusion);i++){if(exclusion[i] in k){b=1}}; if(!(a==1 && b==1)){print}}' ../stringtie_pipeline/correspondingList stringtie.count_5.AS.out > AS_events_withNovelIsoform.tsv

awk '{print $2}' AS_events_withNovelIsoform.tsv|sort | uniq -c

## Get fusion transcripts and their associated reference gene loci
Rscript fusion_tx_detection.R -q stringtie.count_5.strandEnforced.gtf -s ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gtf -t 20 -o fusion_candidates.tsv

## Generate TPM and rawCounts of fusion transcripts
## TPM:
awk 'NR==FNR&&FNR>1{a[$1]=$2}NR>FNR&&FNR==1{print $0"\trefGene"}NR>FNR&&FNR>1{if($1 in a){print $0"\t"a[$1]}}' fusion_candidates.tsv ../expressions/Nanopore/salmonminimap2/Nano_TPM.tsv > fusion_candidates.TPM.tsv
## Raw counts (not from salmon output)
awk 'NR==FNR&&FNR>1{a[$1]=$2}NR>FNR&&FNR==1{print $0"\trefGene"}NR>FNR&&FNR>1{if($1 in a){print $0"\t"a[$1]}}' fusion_candidates.tsv ../expressions/Nanopore/CPM/readsCount_longFormat.out > fusion_candidates.rawCounts.Long.tsv

Rscript long2wide.R
