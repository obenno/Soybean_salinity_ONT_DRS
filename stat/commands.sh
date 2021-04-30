#! /usr/bin/env bash

## Merge minimap2 vs genome bam files into one
samtools merge -@ 20 nanoReads_vs_genome.minimap2.bam ../minimap2_mapping/C08*.bam

## Get identity from bam, reads length from fastq
## BLAST identity was used, perl code is from Heng Li's blog:
## https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
zcat ../fastq/*.fq.gz | awk 'NR==FNR && FNR%4==1{id=substr($1,2); getline; len[id]=length($1)}NR>FNR{print $1"\t"$2"\t"len[$1]}' - <(samtools view -F 0x4 -F 0x100 -F 0x800 nanoReads_vs_genome.minimap2.bam | perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print($F[0], ,"\t", ($l-$n)/$l, "\n")}') > nanoReads_identity_length.tsv

## Plot heat scatter plot with nanoReads_identity_length.tsv
Rscript plot_heatScatter_nanoReads_identity.R nanoReads_identity_length.tsv nanoReads_identity_length.pdf

## Assign reads to isoform according to minimap2 result (against transcriptome)
## Merge minimap2 bam files into one
samtools merge -@ 20 nanoReads_vs_StringTieTranscriptome.minimap2.bam ../minimap2_mapping/NanoReads_vs_StringTieTranscripts.*.sorted.bam
## Generate read length and expected length input file
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_StringTieTranscriptome.minimap2.bam| awk 'NR==FNR{a[$1]=$2}NR>FNR{if($3 in a){print $1"\t"$3}}' ../stringtie_pipeline/correspondingList - | awk 'NR==FNR&& FNR%4==1{id=substr($1,2); getline; len[id]=length($1)}NR>FNR{print $1"\t"$2"\t"len[$1]}' <(zcat ../fastq/*.fq.gz) - | awk 'NR==FNR{a[$1]=length($2)}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2]}' <(fasta_formatter -t -i ../stringtie_pipeline/stringtie.annotated.transcript.fa) - > nanoReads_length_expectedLen.tsv
## Plot heat scatter plot with nanoReads_length_expectedLen.tsv
Rscript plot_heatScatter_nanoReads_length.R nanoReads_length_expectedLen.tsv nanoReads_length_expectedLen.pdf

## Merge vs reference transcriptome bam files into one
samtools merge -@ 20 nanoReads_vs_referenceTranscriptome.minimap2.bam ../minimap2_mapping/NanoReads_vs_Gmax_a2v1_Transcripts.C08*.minimap2.sorted.bam
## Calculate transcript coverage per alignment
## StringTie assembled transcripts
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_StringTieTranscriptome.minimap2.bam | awk 'NR==FNR{a[$1]=length($NF)}NR>FNR{print a[$3]"\t"$0}' <(fasta_formatter -t -i ../stringtie_pipeline/stringtie.annotated.transcript.fa) - | awk 'NR==FNR{a[$1]=$2}NR>FNR{if($4 in a){print "Annotated\t"$0}else{print "Novel\t"$0}}' ../stringtie_pipeline/correspondingList - | perl -ane '$l=0;$l+=$1 while($F[7]=~/(\d+)[MD]/g);print($F[4],"\t", $F[1],"\t",$l/$F[1], "\t", $F[0],"\n")'| sort -k 1,1 | bedtools groupby -g 1 -c 2,3,4 -o distinct,mean,distinct > stringtieTranscriptCoverage.tsv
## Add length groups
awk '$2<=500{print $0"\t(0,500]"}$2>500&&$2<=1000{print $0"\t(500,1000]"}$2>1000&&$2<=1500{print $0"\t(1000,1500]"}$2>1500&&$2<=2000{print $0"\t(1500,2000]"}$2>2000{print $0"\t(2000,max]"}' stringtieTranscriptCoverage.tsv > combined_transcripts_coverage.tsv
## Reference transcripts
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_referenceTranscriptome.minimap2.bam | awk 'NR==FNR{a[$1]=length($NF)}NR>FNR{print a[$3]"\t"$0}' <(fasta_formatter -t -i /home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.transcript.fa) - | perl -ane '$l=0;$l+=$1 while($F[6]=~/(\d+)[MD]/g);print($F[3],"\t", $F[0],"\t",$l/$F[0], "\n")' | sort -k 1,1 | bedtools groupby -g 1 -c 2,3 -o distinct,mean | awk '{print $0"\tReference"}' > referenceTranscriptsCoverage.tsv
## Add groups
awk '$2<=500{print $0"\t(0,500]"}$2>500&&$2<=1000{print $0"\t(500,1000]"}$2>1000&&$2<=1500{print $0"\t(1000,1500]"}$2>1500&&$2<=2000{print $0"\t(1500,2000]"}$2>2000{print $0"\t(2000,max]"}' referenceTranscriptsCoverage.tsv >> combined_transcripts_coverage.tsv
## Plot coverage
Rscript plot_transCoverage.R combined_transcripts_coverage.tsv transcripts_coverage.pdf

####################################

## Get transcript counts for all genes
## StringTie transcriptome
awk '$3=="transcript"{print $12"\t"$10}' ../stringtie_pipeline/stringtie.annotated.gtf | sed 's/\"//g;s/;//g' | awk 'NR==FNR{a[$1]=$2}NR>FNR{if($2 in a){print $0"\tAnnotated"}else{print $0"\tNovel"}}' ../stringtie_pipeline/correspondingList - | sort -k 1,1 -k 2,2 | bedtools groupby -g 1 -c 2,3 -o count,distinct| awk '{split($3, tmp, ","); print $1"\t"$2"\t"tmp[1]}' > gene_transcriptCount.tsv
## Reference transcriptome
awk '$3=="transcript"{print $12"\t"$10}' ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.gene_exons.gtf | sed 's/\"//g; s/;//g' | sort -k 1,1 | bedtools groupby -g 1 -c 2 -o count | awk '{print $1"\t"$2"\tReference"}' >> gene_transcriptCount.tsv
## Plot gene transcript count
Rscript plot_geneTranscriptCounts.R gene_transcriptCount.tsv gene_transcriptCount.pdf
