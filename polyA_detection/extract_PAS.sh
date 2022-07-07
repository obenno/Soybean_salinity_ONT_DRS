#! /usr/bin/env bash


## Extrac PAS site from nanopolish polya tsv ouput
inputTSV=$1
inputBAM=$2
outputFile=$3

## Extract genomic coordinates from bam file
## Filter PAS reads by soft clip: reads with > 10nt soft clip at 3' end were discarded
## Only primary mapping results were used (-F 0x100), SUPPLEMENTARY alignment were also excluded (-F 0x800)
awk '$NF=="PASS"{print $1}' $inputTSV |
    awk 'NR==FNR{a[$1]}NR>FNR{if($1 in a){print}}' - <(samtools view -F 0x4 -F 0x100 -F 0x800 $inputBAM) |
    perl -ane 'if(($F[1] & 0x10) ==0){$F[5]!~/(\d+)S$/; if($1<=10){print join("\t", @F), "\n"}}else{$F[5]!~/^(\d+)S/; if($1<=10){print join("\t", @F), "\n"}}' |
    perl -ane 'if(($F[1] & 0x10)==0){$strand="+";$l=0;$l+=$1 while($F[5]=~/(\d+)[MDN]/g);$pos=$F[3]+$l-1}else{$strand="-";$pos=$F[3]}; print($F[0],"\t",$F[2],"\t",$pos,"\t",$strand, "\n")' |
    awk 'BEGIN{print "ReadID\tPolyA_len\tChr\tPAS_positon\tStrand"}NR==FNR{a[$1]=$9}NR>FNR{print $1"\t"a[$1]"\t"$2"\t"$3"\t"$4}' $inputTSV - > $outputFile

