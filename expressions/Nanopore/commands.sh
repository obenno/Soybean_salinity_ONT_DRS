#! /usr/bin/env bash

## First randomize minimap2 mapping result
bash randomize_bam.sh ../../../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08leaf_0h.minimap2.bam NanoReads_vs_StringTieTranscripts.C08leaf_0h.minimap2.randomized.bam
bash randomize_bam.sh ../../../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08leaf_1h.minimap2.bam NanoReads_vs_StringTieTranscripts.C08leaf_1h.minimap2.randomized.bam
bash randomize_bam.sh ../../../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08root_0h.minimap2.bam NanoReads_vs_StringTieTranscripts.C08root_0h.minimap2.randomized.bam
bash randomize_bam.sh ../../../minimap2_mapping/NanoReads_vs_StringTieTranscripts.C08root_1h.minimap2.bam NanoReads_vs_StringTieTranscripts.C08root_1h.minimap2.randomized.bam

## Salmon alignment based commands
salmon quant -t ../../illumina/stringtie.annotated.transcript.fa -l A --noErrorModel -a NanoReads_vs_StringTieTranscripts.C08leaf_0h.minimap2.randomized.bam -o C08leaf_0h -p 30
salmon quant -t ../../illumina/stringtie.annotated.transcript.fa -l A --noErrorModel -a NanoReads_vs_StringTieTranscripts.C08leaf_1h.minimap2.randomized.bam -o C08leaf_1h -p 30
salmon quant -t ../../illumina/stringtie.annotated.transcript.fa -l A --noErrorModel -a NanoReads_vs_StringTieTranscripts.C08root_0h.minimap2.randomized.bam -o C08root_0h -p 30
salmon quant -t ../../illumina/stringtie.annotated.transcript.fa -l A --noErrorModel -a NanoReads_vs_StringTieTranscripts.C08root_1h.minimap2.randomized.bam -o C08root_1h -p 30

## Plot cor for nano TPM and illumina TPM
## Nano TPM table will also be generated
Rscript plot_TPM_cor.R

## Add ref Gene ID
awk 'NR==FNR{a[$1]=$2}NR>FNR&&FNR>1{if($2 in a){print $1"\t"$2"\t"a[$2]}else{print $1"\t"$2"\tNA"}}' ../../../DASG_PSI/known_genes.lst tx2gene.tsv | awk 'NR==FNR{a[$1]=$3"\t"$2}NR>FNR&&FNR==1{print "RefGene\tGeneID\t"$0}NR>FNR&&FNR>1{print a[$1]"\t"$0}' - Nano_TPM.tsv > Nano_TPM.withRefGeneID.tsv
