#! /usr/bin/env bash

## Activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate AS_env

## These commands were used to get PSI result from AS events
## other than aPAS

## Detect and extract AS events from stringtie gtf
awk 'BEGIN{FS="\t";OFS="\t"}{if($7=="."){$7="+"};print}' ../stringtie/stringtie.count_5.gtf > stringtie.count_5.strandEnforced.gtf
python extract_AS_events.py -i stringtie.count_5.strandEnforced.gtf -o stringtie.count_5.AS.out

## Ensure at least one inclusion isoform and one exclusion isoform has max TPM >= 1
awk 'NR==FNR{a[$1]}NR>FNR{split($7, inclusion, ","); split($8, exclusion, ","); m=0; n=0; for(i=1;i<=length(inclusion);i++){if(inclusion[i] in a){m=1}}; for(i=1;i<=length(exclusion);i++){if(exclusion[i] in a){n=1}}; if(m==1 && n==1){print}}' ../expressions/Nano_maxTPM1.lst stringtie.count_5.AS.out > stringtie.count_5.maxTPM_1.AS.out

## Add chr and strand information to the AS table
awk '$3=="transcript"{print $1"\t"$7"\t"$12}' ../stringtie/stringtie.count_5.gtf | sed 's/\"//g; s/;//'|sort -u |
    awk 'NR==FNR{a[$3]=$1"\t"$2}NR>FNR{print $1"\t"a[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' - stringtie.count_5.maxTPM_1.AS.out > stringtie.count_5.maxTPM_1.AS.withChrStrand.out


## Format AS events inclusion and exclusion count table and calculate PSI value
Rscript format_AS_expr.R

## User rmats to test DS events
conda activate rmats-stat

python rMATS_unpaired.py rmats_input.Leaf_1h_vs_Leaf_0h.tsv rmats_output 40 0.1
mv rmats_output/rMATS_Result_P.txt rmats_output/rMATS_Result_P.Leaf_1h_vs_Leaf_0h.txt
python FDR.py rmats_output/rMATS_Result_P.Leaf_1h_vs_Leaf_0h.txt rmats_output/rMATS_Result_FDR.Leaf_1h_vs_Leaf_0h.txt

python rMATS_unpaired.py rmats_input.Root_1h_vs_Root_0h.tsv rmats_output 40 0.1
mv rmats_output/rMATS_Result_P.txt rmats_output/rMATS_Result_P.Root_1h_vs_Root_0h.txt
python FDR.py rmats_output/rMATS_Result_P.Root_1h_vs_Root_0h.txt rmats_output/rMATS_Result_FDR.Root_1h_vs_Root_0h.txt

python rMATS_unpaired.py rmats_input.Root_0h_vs_Leaf_0h.tsv rmats_output 40 0.1
mv rmats_output/rMATS_Result_P.txt rmats_output/rMATS_Result_P.Root_0h_vs_Leaf_0h.txt
python FDR.py rmats_output/rMATS_Result_P.Root_0h_vs_Leaf_0h.txt rmats_output/rMATS_Result_FDR.Root_0h_vs_Leaf_0h.txt

python rMATS_unpaired.py rmats_input.Root_1h_vs_Leaf_1h.tsv rmats_output 40 0.1
mv rmats_output/rMATS_Result_P.txt rmats_output/rMATS_Result_P.Root_1h_vs_Leaf_1h.txt
python FDR.py rmats_output/rMATS_Result_P.Root_1h_vs_Leaf_1h.txt rmats_output/rMATS_Result_FDR.Root_1h_vs_Leaf_1h.txt

conda deactivate

## Only include all AS events, no PASs
cat AS_event_PSI.tsv > all_events_psi.tsv

## aggregate psi and FDR tables
Rscript aggregate_psi_FDR.R

## calculate events with FDR <= 0.05
for i in all_events.*.FDR.tsv; do count=$(awk '$7<=0.05' $i|wc -l); echo $i $count; done

## Plot DEG and DASG Venn Diagram
awk 'FNR>1{print $1}' ../expressions/*_1hvs0h.gene_up.out ../expressions/*_1hvs0h.gene_down.out | sort -u > DEG.lst
awk 'FNR>1 && $7<=0.05' all_events.Leaf_1h_vs_Leaf_0h.FDR.tsv all_events.Root_1h_vs_Root_0h.FDR.tsv | cut -f 2 | sort -u > DASG.lst

Rscript vennPlot.R -i DEG.lst,DASG.lst -n DEG,DASG -o DEG_vs_DASG.venn.pdf


## Count how many events containing one isofrom completely from novel transcripts
awk 'NR==FNR{k[$1]=$2}NR>FNR&&FNR>1{split($7, inclusion, ","); split($8, exclusion, ","); a=0; b=0; for(i=1;i<=length(inclusion);i++){if(inclusion[i] in k){a=1}}; for(i=1;i<=length(exclusion);i++){if(exclusion[i] in k){b=1}}; if(!(a==1 && b==1)){print}}' ../stringtie/correspondingList stringtie.count_5.AS.out > AS_events_withNovelIsoform.tsv

awk '{print $2}' AS_events_withNovelIsoform.tsv|sort | uniq -c

## Get fusion transcripts and their associated reference gene loci
Rscript fusion_tx_detection.R -q stringtie.count_5.strandEnforced.gtf -s ~/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 -t 20 -o fusion_candidates.tsv

## Generate TPM of fusion transcripts
awk 'NR==FNR&&FNR>1{a[$1]=$2}NR>FNR&&FNR==1{print $0"\trefGene"}NR>FNR&&FNR>1{if($1 in a){print $0"\t"a[$1]}}' fusion_candidates.tsv ../expressions/Nano_TPM.tsv > fusion_candidates.TPM.tsv
