#! /usr/bin/env bash

## Activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate salmon

## First randomize minimap2 mapping result
## then mapping with salmon
for i in ../stringtie/*.vsStringTieTranscriptome.minimap2.bam
do
    sampleBAM=$(basename $i)
    sampleName=${sampleBAM%%.vsStringTieTranscriptome.minimap2.bam}
    nameSortedBAM=$sampleName".vsStringTieTranscriptome.nameSorted.bam"
    samtools sort -@ 40 -n -o $nameSortedBAM $i
    salmon quant -t ../stringtie/stringtie.annotated.transcript.fa -l A --noErrorModel -p 30 -a $nameSortedBAM -o $sampleName
done

## Generate tx2gene.tsv file
awk -F"\t" 'BEGIN{print "TXNAME\tGENEID\tRef_transcript"}$3=="transcript"{split($9,tmp,";"); for(i=1;i<=length(tmp);i++){sub(/^ /, "", tmp[i]); split(tmp[i], k, " "); gsub(/"/,"", k[2]);a[k[1]]=k[2]}; if(a["class_code"]=="="){print a["transcript_id"]"\t"a["gene_id"]"\t"a["cmp_ref"]}else{print a["transcript_id"]"\t"a["gene_id"]"\tNA"}}' ../stringtie/stringtie.annotated.strandCorrected.gtf > tx2gene.tsv

## Generate gene location file
awk '$3=="transcript"{print $1"\t"$4"\t"$5"\t"$7"\t"$12}' ../transdecoder/stringtie.count_5.withCDS.gtf|
    sed 's/\"//g; s/;//' |
    sort -k 5,5 -k 1,1 -k 2,2 -k 3,3 |
    bedtools groupby -g 5 -c 1,2,3,4 -o distinct,min,max,distinct |
    awk '{print $1"\t"$2":"$3"-"$4"("$5")"}' > geneLOC.tsv

## Calculate DEGs and expression table and DTE result
Rscript calculateDEG.R

## Add ref Gene ID
awk 'NR==FNR{a[$1]=$2}NR>FNR&&FNR>1{if($2 in a){print $1"\t"$2"\t"a[$2]}else{print $1"\t"$2"\tNA"}}' ../stringtie/known_genes.lst tx2gene.tsv | awk 'NR==FNR{a[$1]=$3"\t"$2}NR>FNR&&FNR==1{print "RefGene\tGeneID\t"$0}NR>FNR&&FNR>1{print a[$1]"\t"$0}' - Nano_TPM.tsv > Nano_TPM.withRefGeneID.tsv

