#! /usr/bin/env bash

## Commands for m6A sites and ratio stats

## Generate whole genome m6A sites (without depth filtration)
sort -k 3,3 -k 1,1 -k 2,2n -k 4,4 ../C08root_0h/result_final/genome_abandance.0.5.bed | cut -f 1,2,6 | sort -u | sort -k 1,1 -k 2,2n > m6A_wholeGenomeSites.C08root_0h.tsv
sort -k 3,3 -k 1,1 -k 2,2n -k 4,4 ../C08root_1h/result_final/genome_abandance.0.5.bed | cut -f 1,2,6 | sort -u | sort -k 1,1 -k 2,2n > m6A_wholeGenomeSites.C08root_1h.tsv
sort -k 3,3 -k 1,1 -k 2,2n -k 4,4 ../C08leaf_0h/result_final/genome_abandance.0.5.bed | cut -f 1,2,6 | sort -u | sort -k 1,1 -k 2,2n > m6A_wholeGenomeSites.C08leaf_0h.tsv
sort -k 3,3 -k 1,1 -k 2,2n -k 4,4 ../C08leaf_1h/result_final/genome_abandance.0.5.bed | cut -f 1,2,6 | sort -u | sort -k 1,1 -k 2,2n > m6A_wholeGenomeSites.C08leaf_1h.tsv
cat m6A_wholeGenomeSites.*.tsv | sort -u | sort -k 1,1 -k 2,2n > m6A_wholeGenomeSites.tsv

## Plot m6A sites ratio distribution, m6A sites were also extracted (depth >=20)
awk '{split($1, tmp, "|"); gene=tmp[1]; chr=tmp[2]; for(i=2;i<=NF;i++){split($i, d, "|"); pos=d[1]; methyl_count=d[2]; total_count=d[3]; ratio=d[4]; print gene"\t"chr"\t"pos"\t"ratio"\t"methyl_count"\t"total_count}}' ../C08leaf_0h/result_final/ratio.0.5.tsv | awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2"\t"$3]"\t"$4"\t"$5"\t"$6}' m6A_wholeGenomeSites.tsv - > m6A_ratio.C08leaf_0h.tsv
awk '{split($1, tmp, "|"); gene=tmp[1]; chr=tmp[2]; for(i=2;i<=NF;i++){split($i, d, "|"); pos=d[1]; methyl_count=d[2]; total_count=d[3]; ratio=d[4]; print gene"\t"chr"\t"pos"\t"ratio"\t"methyl_count"\t"total_count}}' ../C08leaf_1h/result_final/ratio.0.5.tsv | awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2"\t"$3]"\t"$4"\t"$5"\t"$6}' m6A_wholeGenomeSites.tsv - > m6A_ratio.C08leaf_1h.tsv
awk '{split($1, tmp, "|"); gene=tmp[1]; chr=tmp[2]; for(i=2;i<=NF;i++){split($i, d, "|"); pos=d[1]; methyl_count=d[2]; total_count=d[3]; ratio=d[4]; print gene"\t"chr"\t"pos"\t"ratio"\t"methyl_count"\t"total_count}}' ../C08root_0h/result_final/ratio.0.5.tsv | awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2"\t"$3]"\t"$4"\t"$5"\t"$6}' m6A_wholeGenomeSites.tsv - > m6A_ratio.C08root_0h.tsv
awk '{split($1, tmp, "|"); gene=tmp[1]; chr=tmp[2]; for(i=2;i<=NF;i++){split($i, d, "|"); pos=d[1]; methyl_count=d[2]; total_count=d[3]; ratio=d[4]; print gene"\t"chr"\t"pos"\t"ratio"\t"methyl_count"\t"total_count}}' ../C08root_1h/result_final/ratio.0.5.tsv | awk 'NR==FNR{a[$1"\t"$2]=$3}NR>FNR{print $1"\t"$2"\t"$3"\t"a[$2"\t"$3]"\t"$4"\t"$5"\t"$6}' m6A_wholeGenomeSites.tsv - > m6A_ratio.C08root_1h.tsv
## Plot bar plots for calculating m6A sites
Rscript plot_m6A_sites_bar.R
## Plot venn diagram of m6A sites for four samples
for i in m6A_ratio.*.tsv; do
    awk '{print $1":"$2":"$3}' $i > ${i%%.tsv}".lst"
done

Rscript vennPlot.R -i m6A_ratio.C08root_0h.lst,m6A_ratio.C08root_1h.lst,m6A_ratio.C08leaf_0h.lst,m6A_ratio.C08leaf_1h.lst \
        -n C08root_0h,C08root_1h,C08leaf_0h,C08leaf_1h -o m6A_sites.venn.pdf \
        --width 2.5 --height 2.5
rm VennDiagram*.log

## Plot ratio level as boxplot
Rscript plot_m6A_ratio_boxplot.R
## Generate gene TPM table
awk 'NR==FNR{a[$1]=$2}NR>FNR && FNR>1{print a[$1]"\t"$0}' ../../polyA_detection/trans_gene.tsv ../../expressions/Nanopore/salmonminimap2/Nano_TPM.long.tsv > Nanopore_salmonminimap2_TPM.long.tsv
## Plot ratio vs gene expression
Rscript plot_m6A_ratio_vs_geneExpr.R
## Detect differential methylated sites
## Fisher exact test were employed
## scatter plot for all common sites were generated
Rscript detect_differential_methylated_sites.R

## Preprare up/down m6A gene list for GO enrichment
awk '$9-$8>=0.1' C08leaf_dm6A.out.tsv | awk '{split($1, tmp, ":"); print tmp[1]}'| sort -u | grep -f - ../../DASG_PSI/known_genes.lst |awk '{print $2}'| sort -u > C08leaf_dm6A.up.refGene.lst
awk '$9-$8<=-0.1' C08leaf_dm6A.out.tsv | awk '{split($1, tmp, ":"); print tmp[1]}'| sort -u | grep -f - ../../DASG_PSI/known_genes.lst |awk '{print $2}'| sort -u > C08leaf_dm6A.down.refGene.lst
awk '$9-$8>=0.1' C08root_dm6A.out.tsv | awk '{split($1, tmp, ":"); print tmp[1]}'| sort -u | grep -f - ../../DASG_PSI/known_genes.lst |awk '{print $2}'| sort -u > C08root_dm6A.up.refGene.lst
awk '$9-$8<=-0.1' C08root_dm6A.out.tsv | awk '{split($1, tmp, ":"); print tmp[1]}'| sort -u | grep -f - ../../DASG_PSI/known_genes.lst |awk '{print $2}'| sort -u > C08root_dm6A.down.refGene.lst

## Generate m6A sites bed files
awk '{print $2"\t"$3-1"\t"$3}' m6A_ratio.C08leaf_0h.tsv > m6A_sites.C08leaf_0h.bed
awk '{print $2"\t"$3-1"\t"$3}' m6A_ratio.C08leaf_1h.tsv > m6A_sites.C08leaf_1h.bed
awk '{print $2"\t"$3-1"\t"$3}' m6A_ratio.C08root_0h.tsv > m6A_sites.C08root_0h.bed
awk '{print $2"\t"$3-1"\t"$3}' m6A_ratio.C08root_1h.tsv > m6A_sites.C08root_1h.bed

## Select differential m6A site and check their distribution in gene regions
## Leaf:
awk 'NR>1 && $9!="NA" && $8!="NA" && $9-$8>=0.1' C08leaf_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08leaf_m6A_upSites.0.1.bed
awk 'NR>1 && $9!="NA" && $8!="NA" && $9-$8<=-0.1' C08leaf_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08leaf_m6A_downSites.0.1.bed
awk 'NR>1 && $9!="NA" && $8!="NA" && ($9-$8<0.1||$9-$8>-0.1)' C08leaf_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08leaf_m6A_unchangedSites.0.1.bed

Rscript plot_m6A_geneFeatures_density.R -i ../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf -s C08leaf_m6A_upSites.0.1.bed,C08leaf_m6A_unchangedSites.0.1.bed,C08leaf_m6A_downSites.0.1.bed -n C08leaf_up,C08leaf_unchanged,C08leaf_down -t 30 -o C08leaf_dm6A_distribution.delta0.1.withUnchanged.pdf
## Root:
awk 'NR>1 && $9!="NA" && $8!="NA" && $9-$8>=0.1' C08root_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08root_m6A_upSites.0.1.bed
awk 'NR>1 && $9!="NA" && $8!="NA" && $9-$8<=-0.1' C08root_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08root_m6A_downSites.0.1.bed
awk 'NR>1 && $9!="NA" && $8!="NA" && ($9-$8<0.1||$9-$8>-0.1)' C08root_dm6A.out.tsv | awk '{split($1,tmp, ":"); print tmp[2]"\t"tmp[3]-1"\t"tmp[3]}'| sort -u | sort -k 1,1 -k 2,2n > C08root_m6A_unchangedSites.0.1.bed

Rscript plot_m6A_geneFeatures_density.R -i ../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf -s C08root_m6A_upSites.0.1.bed,C08root_m6A_unchangedSites.0.1.bed,C08root_m6A_downSites.0.1.bed -n C08root_up,C08root_unchanged,C08root_down -t 30 -o C08root_dm6A_distribution.delta0.1.withUnchanged.pdf

## Plot methylation ratio boxplot of up, down, unchanged sites
## grouped by gene features (5' UTR, CDS, 3' UTR)
Rscript plot_m6A_ratio_geneFeatures.boxplot.R

