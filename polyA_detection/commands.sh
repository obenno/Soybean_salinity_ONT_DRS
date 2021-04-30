#! /usr/bin/env bash

## Commands for PAS detection and polyA length estimation

nanopolish polya --threads=8 --reads=../fastq/C08leaf_0h.fq.gz --bam=../minimap2_mapping/C08leaf_0h.minimap2.bam --genome=/home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa > C08leaf_0h.polya.tsv

nanopolish polya --threads=8 --reads=../fastq/C08leaf_1h.fq.gz --bam=../minimap2_mapping/C08leaf_1h.minimap2.bam --genome=/home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa > C08leaf_1h.polya.tsv

nanopolish polya --threads=8 --reads=../fastq/C08root_0h.fq.gz --bam=../minimap2_mapping/C08root_0h.minimap2.bam --genome=/home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa > C08root_0h.polya.tsv

nanopolish polya --threads=8 --reads=../fastq/C08root_1h.fq.gz --bam=../minimap2_mapping/C08root_1h.minimap2.bam --genome=/home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa > C08root_1h.polya.tsv

bash extract_PAS.sh C08leaf_0h.polya.tsv ../minimap2_mapping/C08leaf_0h.minimap2.bam C08leaf_0h.PAS.out
bash extract_PAS.sh C08leaf_1h.polya.tsv ../minimap2_mapping/C08leaf_1h.minimap2.bam C08leaf_1h.PAS.out
bash extract_PAS.sh C08root_0h.polya.tsv ../minimap2_mapping/C08root_0h.minimap2.bam C08root_0h.PAS.out
bash extract_PAS.sh C08root_1h.polya.tsv ../minimap2_mapping/C08root_1h.minimap2.bam C08root_1h.PAS.out

awk 'FNR>1{print}' C08leaf_0h.PAS.out C08leaf_1h.PAS.out C08root_0h.PAS.out C08root_1h.PAS.out > combined.PAS.out
awk '$3=="transcript"{print $1"\t"$4"\t"$5"\t"$10"\t"$12}' ../stringtie_pipeline/stringtie.annotated.gtf |sed 's/\"//g; s/;//g' | sort -k 5,5 -k 1,1 -k 2,2 -k 3,3 | bedtools groupby -g 5 -c 1,2,3,4 -o distinct,min,max,distinct | awk '{split($5,tmp,","); for(i=1;i<=length(tmp);i++){print tmp[i]"\t"$1"\t"$2"\t"$3"\t"$4}}' > trans_gene.tsv

## PASs will be merged if their distance is less than 24 nt
python tidy_PAS.py -b ../stat/nanoReads_vs_StringTieTranscriptome.minimap2.bam -p combined.PAS.out -g trans_gene.tsv -d 24 -o collapsed_gene_PAS_reads.tsv

## Filter PAS with transcript read counts (5)
awk 'NR==FNR{a[$1]}NR>FNR{f=0; split($5, trans, ","); for(i=1;i<=length(trans);i++){if(trans[i] in a){f=1}}; if(f==1){print}}' ../stat/stringtie.transcript.counts_more_than_five.lst collapsed_gene_PAS_reads.tsv > collapsed_gene_PAS_reads.transCount5.tsv

## Filter PAS with PAS read counts(5)
awk '{n=0; split($6, read_group, ","); for(i=1;i<=length(read_group);i++){split(read_group[i], tmp, "|"); n+=length(tmp)}; if(n>=5){print}}' collapsed_gene_PAS_reads.transCount5.tsv > collapsed_gene_PAS_reads.transCount5.PASreads5.tsv

## Generate counts table of PASs
awk 'FILENAME==ARGV[1]{lib_A[$1]}FILENAME==ARGV[2]{lib_B[$1]}FILENAME==ARGV[3]{lib_C[$1]}FILENAME==ARGV[4]{lib_D[$1]}FILENAME==ARGV[5]{a=0;b=0;c=0;d=0;split($6, read_group, ","); for(i=1;i<=length(read_group);i++){split(read_group[i], tmp, "|"); for(k=1;k<=length(tmp);k++){if(tmp[k] in lib_A){a++}else if(tmp[k] in lib_B){b++}else if(tmp[k] in lib_C){c++}else if(tmp[k] in lib_D){d++}}}; print $1"\t"$2"\t"$3"\t"$4"\t"a"\t"b"\t"c"\t"d}' C08leaf_0h.PAS.out C08leaf_1h.PAS.out C08root_0h.PAS.out C08root_1h.PAS.out collapsed_gene_PAS_reads.transCount5.PASreads5.tsv > collapsed_gene_PAS_reads.countsTable.tsv

## Add strand infor to PAS count table
awk '$3=="transcript"{print $12"\t"$7}' ../stringtie_pipeline/stringtie.annotated.gtf | sed 's/\"//g; s/;//g' | sort -u | awk 'NR==FNR{a[$1]=$2}NR>FNR{print $1"\t"$2"\t"$3"\t"$4"\t"a[$1]"\t"$5"\t"$6"\t"$7"\t"$8}' - collapsed_gene_PAS_reads.countsTable.tsv > collapsed_gene_PAS_reads.strand.countsTable.tsv

## Select the most abundant two PASs and calculate PSI
Rscript select_PAS.R

## Generate PAS bed file, the selected distal and proximal
## PASs for gene with multiple PASs will be colored by dark purple
awk 'NR==FNR && $2>1{a[$3"\t"$4"\t"$5]}NR>FNR{n=$6+$7+$8+$9; if($2"\t"$3"\t"$4 in a){color="88,24,69"}else{color="."};print $2"\t"$3-1"\t"$4"\t"n"\t"0"\t"$5"\t"$3-1"\t"$4"\t"color"\t1\t"$4-$3+1",\t0,"}' selected_PAS_cpm.tsv collapsed_gene_PAS_reads.strand.countsTable.tsv > filtered_PAS.bed

## Format PASs' CPM result
cut -f 1,2,3-6,11-14 selected_PAS_cpm.tsv | awk 'NR>1 && $2 >1'| sort -k 1,1 -k 4,4n | awk '{print $1"\t"$3":"$4"-"$5":"$6"\t"$7"|"$8"|"$9"|"$10}'| bedtools groupby -g 1 -c 2,3 -o collapse | awk '{split($2, exon, ","); split($3, expr, ","); strand=substr(exon[1],length(exon[1]));if(strand=="+"){print $1"\taPAS\t"exon[2]"\t"exon[1]"\t"expr[2]"\t"expr[1]}else{print $1"\taPAS\t"exon[1]"\t"exon[2]"\t"expr[1]"\t"expr[2]}}' > aPAS_events_cpm.tsv

## Generate transcript polyA length profile
bash extract_polyA_length.sh C08leaf_0h.PAS.out ../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08leaf_0h.minimap2.bam transcript_polyA_length.C08leaf_0h.tsv
bash extract_polyA_length.sh C08leaf_1h.PAS.out ../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08leaf_1h.minimap2.bam transcript_polyA_length.C08leaf_1h.tsv
bash extract_polyA_length.sh C08root_0h.PAS.out ../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08root_0h.minimap2.bam transcript_polyA_length.C08root_0h.tsv
bash extract_polyA_length.sh C08root_1h.PAS.out ../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08root_1h.minimap2.bam transcript_polyA_length.C08root_1h.tsv

## Generate violin plot for polyA length of each sample
## and perform Mann–Whitney U test
Rscript process_polyA_length.R

## Compare detected PAS with stringtie transcripts end and annotated transcript end
awk '$3=="transcript"{if($7=="+"){print $1"\t"$5"\t"$5+1"\t"$10"_end\t.\t"$7}else{print $1"\t"$4-1"\t"$4"\t"$10"_end\t.\t"$7}}' ../alternative_splicing/stringtie.count_5.gtf | sed 's/\"//g; s/;//g' | sort -k 1,1 -k 2,2n | bedtools merge -s -c 4,5,6 -o distinct| sort -k 1,1 -k 2,2n > stringtie.count_5.three_end.bed

awk '$3=="transcript"{if($7=="+"){print $1"\t"$5"\t"$5+1"\t"$10"_end\t.\t"$7}else{print $1"\t"$4-1"\t"$4"\t"$10"_end\t.\t"$7}}' ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gtf |sed 's/\"//g; s/;//g' | sort -k 1,1 -k 2,2n | bedtools merge -s -c 4,5,6 -o distinct| sort -k 1,1 -k 2,2n > Gmax_a2v1.transcripts.three_end.bed

bedtools closest -k 1 -D b -s -a <(cut -f 1-6 filtered_PAS.bed|sort -k 1,1 -k 2,2n ) -b Gmax_a2v1.transcripts.three_end.bed | cut -f 1-6,13 > nanoPAS_vs_Gmax_a2v1_three_end.bed

bedtools closest -k 1 -D b -s -a <(cut -f 1-6 filtered_PAS.bed|sort -k 1,1 -k 2,2n ) -b stringtie.count_5.three_end.bed | cut -f 1-6,13 > nanoPAS_vs_StringTie_count5_three_end.bed

## Plot distribution of distance between identified PASs and annotated transcripts ends
Rscript plot_nanoPAS_vs_AnnotationDistance.R

## Add gene annotation for alternative polyA length transcripts
awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR==1{print "gene\t"$0}NR>FNR && FNR>1{print a[$1]"\t"$0}' trans_gene.tsv alternative_polyA_length.C08leaf.tsv | awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR==1{print "refGene\t"$0}NR>FNR && FNR>1{if($1 in a){print a[$1]"\t"$0}else{print "NA\t"$0}}' ../DASG_PSI/known_genes.lst - | awk -F"\t" 'NR==FNR{split($1, tmp, "."); a[tmp[1]"."tmp[2]]=$3}NR>FNR && FNR==1{print $0"\tAnno"}NR>FNR && FNR>1{if($1 in a){print $0"\t"a[$1]}else{print $0"\tNA"}}' ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.defline.txt - > alternative_polyA_length.C08leaf.withAnno.tsv

awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR==1{print "gene\t"$0}NR>FNR && FNR>1{print a[$1]"\t"$0}' trans_gene.tsv alternative_polyA_length.C08root.tsv | awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR==1{print "refGene\t"$0}NR>FNR && FNR>1{if($1 in a){print a[$1]"\t"$0}else{print "NA\t"$0}}' ../DASG_PSI/known_genes.lst - | awk -F"\t" 'NR==FNR{split($1, tmp, "."); a[tmp[1]"."tmp[2]]=$3}NR>FNR && FNR==1{print $0"\tAnno"}NR>FNR && FNR>1{if($1 in a){print $0"\t"a[$1]}else{print $0"\tNA"}}' ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.defline.txt - > alternative_polyA_length.C08root.withAnno.tsv
## Add transcript expr to the table
awk 'NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5}NR>FNR&&FNR==1{print $0"\tC08leaf_0h_TPM\tC08leaf_1h_TPM\tC08root_0h_TPM\tC08root_1h_TPM"}NR>FNR&&FNR>1{print $0"\t"a[$3]}' ../expressions/Nanopore/salmonminimap2/Nano_TPM.tsv alternative_polyA_length.C08leaf.withAnno.tsv > alternative_polyA_length.C08leaf.withAnno.withExpr.tsv
awk 'NR==FNR{a[$1]=$2"\t"$3"\t"$4"\t"$5}NR>FNR&&FNR==1{print $0"\tC08leaf_0h_TPM\tC08leaf_1h_TPM\tC08root_0h_TPM\tC08root_1h_TPM"}NR>FNR&&FNR>1{print $0"\t"a[$3]}' ../expressions/Nanopore/salmonminimap2/Nano_TPM.tsv alternative_polyA_length.C08root.withAnno.tsv > alternative_polyA_length.C08root.withAnno.withExpr.tsv
