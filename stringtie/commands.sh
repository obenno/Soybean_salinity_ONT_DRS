#! /usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate Nanopore

function runStringtie {
    local inputBAM=$1
    local inputFileName=$(basename $inputBAM)
    local outputGTF=${inputFileName%%.vsGenome.minimap2.bam}".stringtie.gtf"
    stringtie -G /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 \
              -p 40 -L \
              -o stringtie_singleSample/$outputGTF $inputBAM
}

## StringTie commands
for i in ../minimap2_mapping/*.vsGenome.minimap2.primary.bam;
do
    runStringtie $i
done
## StringTie merge
stringtie --merge -G /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 \
          -F 0 -T 0 -i -o merged.gtf stringtie_singleSample/*.gtf
## Collect isoforms not included in stringtie merge
gffcompare -r merged.gtf stringtie_singleSample/*.gtf
## Comapre and annotate assembled transcripts with annotation
gffcompare -r /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 gffcmp.combined.gtf -o stringtie
## Enforce transcripts with unknown strand info as "+"
awk 'BEGIN{FS="\t";OFS="\t"}{if($7!="+" && $7!="-"){$7="+"}; print}' stringtie.annotated.gtf > stringtie.annotated.strandCorrected.gtf
## Generate transcripts fasta file
gffread -g /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa \
        -w stringtie.annotated.transcript.fa stringtie.annotated.strandCorrected.gtf

## Map raw reads to StringTie transcripts
function map2stringtie {
    inputFastq=$1
    sampleName=$(basename $inputFastq)
    outBAM=${sampleName%%.fq.gz}".vsStringTieTranscriptome.minimap2.bam"
    minimap2 -t 30 -ax map-ont -uf --secondary=no stringtie.annotated.transcript.fa \
             $inputFastq |
        samtools view -F 0x4 -F 0x10 -u |
        samtools sort -@ 30 -l 9 > $outBAM
}

for i in ../../fastq_input/*.fq.gz;
do
    map2stringtie $i;
done

## Merge raw reads vs StringTie transcripts bam files
samtools merge -@ 40 -f nanoReads_vs_StringTieTranscriptome.minimap2.bam *.vsStringTieTranscriptome.minimap2.bam

## Filter stringtie transcript by reads count >= 5
samtools view -F 0x4 -F 0x100 -F 0x800 -F 0x10 nanoReads_vs_StringTieTranscriptome.minimap2.bam |
    awk '{print $1"\t"$3}' | sort -k 2,2 | bedtools groupby -g 2 -c 1 -o count |
    awk '$2>=5{print $1"\t"$2}' > stringtie.transcript.counts_more_than_five.lst

awk 'NR==FNR{a["\""$1"\";"]}NR>FNR{if($10 in a){print}}' stringtie.transcript.counts_more_than_five.lst stringtie.annotated.gtf > stringtie.count_5.gtf

## Prepare colored bed file for jbrowse and IGV

## Generate correspondingList for annotated transcripts
## For non filtered transcripts gtf
awk -F"\t" '$3=="transcript"{split($9,tmp,"; "); for(i=1;i<=length(tmp);i++){split(tmp[i],k," "); a[k[1]]=k[2]};if(a["class_code"]=="\"=\""){print a["transcript_id"]"\t"a["cmp_ref"]}}' stringtie.annotated.gtf |  sed 's/\"//g' > correspondingList
## Generate known_genes.lst
awk '$3=="transcript"{print $10"\t"$12}' stringtie.annotated.gtf |
    sed 's/\"//g; s/;//g' | awk 'NR==FNR{a[$1]=$2}NR>FNR{if($1 in a){print $2"\t"$1"\t"a[$1]}}' correspondingList - |
    awk '{split($3, tmp, "."); print $1"\t"tmp[1]"."tmp[2]}'| sort -u > known_genes.lst

## Convert gtf to genePred
gtfToGenePred -genePredExt stringtie.annotated.gtf stringtie.annotated.gp
## Convert genePred to bed20
genePredToBigGenePred stringtie.annotated.gp stdout | sort -k1,1 -k2,2n > stringtie.annotated.bed20
## add color code for bed20
awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1]}NR>FNR{if(!($4 in a)){if($6=="+"){$9="239,138,98"}else{$9="103,169,207"}}; print}' correspondingList stringtie.annotated.bed20 > stringtie.annotated.colored.bed20

## For filtered transcripts gtf
gtfToGenePred -genePredExt stringtie.count_5.gtf stringtie.count_5.gp
genePredToBigGenePred stringtie.count_5.gp stdout | sort -k1,1 -k2,2n > stringtie.count_5.bed20
awk -F"\t" 'BEGIN{OFS="\t"}NR==FNR{a[$1]}NR>FNR{if(!($4 in a)){if($6=="+"){$9="239,138,98"}else{$9="103,169,207"}}; print}' correspondingList stringtie.count_5.bed20 > stringtie.count_5.colored.bed20

## Generate gene coordinates with transcript ID
awk '$3=="transcript"{print $1"\t"$4"\t"$5"\t"$10"\t"$12}' stringtie.annotated.gtf |
    sed 's/\"//g; s/;//g' | sort -k 5,5 -k 1,1 -k 2,2 -k 3,3 |
    bedtools groupby -g 5 -c 1,2,3,4 -o distinct,min,max,distinct |
    awk '{split($5,tmp,","); for(i=1;i<=length(tmp);i++){print tmp[i]"\t"$1"\t"$2"\t"$3"\t"$4}}' > trans_gene.tsv

## Generate gene bed file for nanom6A
awk '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' stringtie.annotated.strandCorrected.gtf |
    sed 's/\"//g; s/;//' | sort -k 4,4 -k 1,1 -k 2,2n -k 3,3n |
    bedtools groupby -g 4 -c 1,2,3,5,6 -o distinct,min,max,distinct,distinct |
    awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6}' > stringtie.annotated.gene.bed