#! /usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate polya

## The fast5 needs to be indexed by fastq
## nanopolish index -d Leaf_0h_R1 ../fastq_input/Leaf_0h_R1.fq.gz

## Commands for PAS detection and polyA length estimation
## bam index should be prepared firstly
function polya {
    local inputBAM=$1
    local inputFastq=$2
    local fileName=$(basename $inputFastq)
    local output=${fileName%%.fq.gz}".polya.tsv"
    local log=${fileName%%.fq.gz}".polya.log"
    nanopolish polya --threads=20 --reads=$inputFastq \
               --bam=$inputBAM \
               --genome=/home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa \
               > $output 2> $log
}

## ensure vbz_compression plugin works
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin

## generate polya commands and use parallel
ls ../../fastq_input/*.fq.gz |
    awk '{split($1, tmp, "/"); fileName=tmp[4]; sample=substr(fileName, 1, index(fileName, ".fq.gz")-1); print "nanopolish polya --threads=5 --reads="$1" --bam=../minimap2_mapping/"sample".vsGenome.minimap2.primary.bam  --genome=/home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa > "sample".polya.tsv 2> "sample".polya.log"}' |
    parallel -P 12

conda activate samtools

for i in *.polya.tsv
do
    bash extract_PAS.sh $i ../minimap2_mapping/${i%%.polya.tsv}".vsGenome.minimap2.primary.bam" ${i%%.polya.tsv}".PAS.out"
done

