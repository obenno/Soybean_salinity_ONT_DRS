#! /usr/bin/env bash

## Calculate raw reads numbers,
## reads mean length, max length, and length N50
function len_stat {
    local inputFastq=$1
    local fileName=$(basename $inputFastq)
    local sampleName=${fileName%%.fq.gz}
    zcat $inputFastq |
        awk '
        BEGIN{len_sum=0; len_max=0; k=0}
        {
            if(NR%4==1){
                getline;
                k+=1
                a[k]=length($1)
                if(length($1)>len_max){len_max=length($1)};
                len_sum+=length($1)
            }
        }
        END{
            num=asort(a);
            s=0
            n50=0
            for(i=1;i<=num;i++){
               s+=a[i]
               if(s>=len_sum/2){
                  n50 = a[i]
                  break
               }
            }
            len_mean=len_sum/num; print "'$sampleName'\t"num"\t"n50"\t"len_mean"\t"len_max}
        '
}

rm -f rawReadLength_stat
for i in ../../fastq_input/*.fq.gz;
do
    len_stat $i >> rawReadLength_stat
done

## Merge minimap2 vs genome bam files into one
samtools merge -@ 20 nanoReads_vs_genome.minimap2.primary.bam \
         ../minimap2_mapping/*.vsGenome.minimap2.primary.bam

## Get identity from bam, reads length from fastq
## BLAST identity was used, perl code is from Heng Li's blog:
## https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
zcat ../../fastq_input/*.fq.gz |
    awk '
    NR==FNR && FNR%4==1{id=substr($1,2); getline; len[id]=length($1)}
    NR>FNR{print $1"\t"$2"\t"len[$1]}
    ' - <(samtools view -F 0x4 -F 0x100 -F 0x800 nanoReads_vs_genome.minimap2.primary.bam |
              perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print($F[0], "\t", ($l-$n)/$l, "\n")}') > nanoReads_identity_length.tsv

## Plot heat scatter plot with nanoReads_identity_length.tsv (require large RAM)
Rscript plot_heatScatter_nanoReads_identity.R nanoReads_identity_length.tsv nanoReads_identity_length.pdf

## Calculate average mapping identity and mapped reads for each sample
rm -f eachSample_vsGenome_mappingIdentity
for i in ../minimap2_mapping/*.vsGenome.minimap2.primary.bam;
do
    fileName=$(basename $i)
    samtools view -@ 20 -F 0x4 -F 0x100 -F 0x800 $i |
        perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print($F[0], "\t", ($l-$n)/$l, "\n")}' |
        awk 'BEGIN{b=0}{a[$1];b+=$2}END{print "'$fileName'"; print "Mapped Reads: "length(a); print "Average Identity: "b/NR;}' >> eachSample_vsGenome_mappingIdentity
done

rm -f eachSample_vsRefTranscriptome_mappingIdentity
for i in ../minimap2_mapping/*.vsTranscriptome.minimap2.bam;
do
    fileName=$(basename $i)
    samtools view -@ 20 -F 0x4 -F 0x100 -F 0x800 -F 0x10 $i |
        perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print($F[0], "\t", ($l-$n)/$l, "\n")}' |
        awk 'BEGIN{b=0}{a[$1];b+=$2}END{print "'$fileName'"; print "Mapped Reads: "length(a); print "Average Identity: "b/NR;}' >> eachSample_vsRefTranscriptome_mappingIdentity
done

rm -f eachSample_vsStringTieTranscriptome_mappingIdentity
for i in ../stringtie/*.vsStringTieTranscriptome.minimap2.bam;
do
    fileName=$(basename $i)
    samtools view -@ 20 -F 0x4 -F 0x100 -F 0x800 -F 0x10 $i |
        perl -ane 'if(/NM:i:(\d+)/){$n=$1;$l=0;$l+=$1 while/(\d+)[MID]/g;print($F[0], "\t", ($l-$n)/$l, "\n")}' |
        awk 'BEGIN{b=0}{a[$1];b+=$2}END{print "'$fileName'"; print "Mapped Reads: "length(a); print "Average Identity: "b/NR;}' >> eachSample_vsStringTieTranscriptome_mappingIdentity
done

## Assign reads to isoform according to minimap2 result (against transcriptome)
## Merged minimap2 bam files were copied to here and indexed
ln -f -s ../stringtie/nanoReads_vs_StringTieTranscriptome.minimap2.bam ./
samtools index nanoReads_vs_StringTieTranscriptome.minimap2.bam

samtools merge -@ 40 -f nanoReads_vs_referenceTranscriptome.minimap2.bam ../minimap2_mapping/*.vsTranscriptome.minimap2.bam
samtools index nanoReads_vs_referenceTranscriptome.minimap2.bam

## Generate read length and expected length input file
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_referenceTranscriptome.minimap2.bam |
    awk 'NR==FNR&& FNR%4==1{id=substr($1,2); getline; len[id]=length($1)}NR>FNR{print $1"\t"$3"\t"len[$1]}' <(zcat ../../fastq_input/*.fq.gz) - |
    awk 'NR==FNR{a[$1]=length($NF)}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2]}' <(fasta_formatter -t -i /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.transcript.fa) - > nanoReads_length_expectedLen.tsv

## Plot heat scatter plot with nanoReads_length_expectedLen.tsv
Rscript plot_heatScatter_nanoReads_length.R nanoReads_length_expectedLen.tsv nanoReads_length_expectedLen.pdf

## Calculate transcript coverage per alignment
## StringTie assembled transcripts
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_StringTieTranscriptome.minimap2.bam |
    awk 'NR==FNR{a[$1]=length($NF)}NR>FNR{print a[$3]"\t"$0}' <(fasta_formatter -t -i ../stringtie/stringtie.annotated.transcript.fa) - |
    awk 'NR==FNR{a[$1]=$2}NR>FNR{if($4 in a){print "Annotated\t"$0}else{print "Novel\t"$0}}' ../stringtie/correspondingList - |
    perl -ane '$l=0;$l+=$1 while($F[7]=~/(\d+)[MD]/g);print($F[4],"\t", $F[1],"\t",$l/$F[1], "\t", $F[0],"\n")'|
    sort -k 1,1 | bedtools groupby -g 1 -c 2,3,4 -o distinct,mean,distinct > stringtieTranscriptCoverage.tsv
## Add length groups
awk '$2<=500{print $0"\t(0,500]"}$2>500&&$2<=1000{print $0"\t(500,1000]"}$2>1000&&$2<=1500{print $0"\t(1000,1500]"}$2>1500&&$2<=2000{print $0"\t(1500,2000]"}$2>2000{print $0"\t(2000,max]"}' stringtieTranscriptCoverage.tsv > combined_transcripts_coverage.tsv
## Reference transcripts
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_referenceTranscriptome.minimap2.bam |
    awk 'NR==FNR{a[$1]=length($NF)}NR>FNR{print a[$3]"\t"$0}' <(fasta_formatter -t -i /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.transcript.fa) - |
    perl -ane '$l=0;$l+=$1 while($F[6]=~/(\d+)[MD]/g);print($F[3],"\t", $F[0],"\t",$l/$F[0], "\n")' |
    sort -k 1,1 | bedtools groupby -g 1 -c 2,3 -o distinct,mean |
    awk '{print $0"\tReference"}' > referenceTranscriptsCoverage.tsv
## Add groups
awk '$2<=500{print $0"\t(0,500]"}$2>500&&$2<=1000{print $0"\t(500,1000]"}$2>1000&&$2<=1500{print $0"\t(1000,1500]"}$2>1500&&$2<=2000{print $0"\t(1500,2000]"}$2>2000{print $0"\t(2000,max]"}' referenceTranscriptsCoverage.tsv >> combined_transcripts_coverage.tsv
## Plot coverage
Rscript plot_transCoverage.R combined_transcripts_coverage.tsv transcripts_coverage.pdf

## Generate transcript coverage plot with filtered transcripts (counts >= 5)
awk 'NR==FNR{a[$1]}NR>FNR{if(($1 in a)||($1 ~/^Glyma/)){print}}' ../stringtie/stringtie.transcript.counts_more_than_five.lst combined_transcripts_coverage.tsv > combined_transcripts_coverage.count5.tsv

Rscript plot_transCoverage.R combined_transcripts_coverage.count5.tsv transcripts_coverage.count5.pdf

####################################

## Get transcript counts for all genes
## StringTie transcriptome
awk '$3=="transcript"{print $12"\t"$10}' ../stringtie/stringtie.annotated.gtf | sed 's/\"//g;s/;//g' |
    awk 'NR==FNR{a[$1]=$2}NR>FNR{if($2 in a){print $0"\tAnnotated"}else{print $0"\tNovel"}}' ../stringtie/correspondingList - |
    sort -k 1,1 -k 2,2 | bedtools groupby -g 1 -c 2,3 -o count,distinct|
    awk '{split($3, tmp, ","); print $1"\t"$2"\t"tmp[1]}' > gene_transcriptCount.tsv
## Reference transcriptome
awk '$3=="transcript"{print $12"\t"$10}' /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gtf |
    sed 's/\"//g; s/;//g' | sort -k 1,1 |
    bedtools groupby -g 1 -c 2 -o count | awk '{print $1"\t"$2"\tReference"}' >> gene_transcriptCount.tsv
## Plot gene transcript count
Rscript plot_geneTranscriptCounts.R gene_transcriptCount.tsv gene_transcriptCount.pdf
