#! /usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate Nanopore

set -euo pipefail

## ensure vbz_compression plugin works
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin

## Add basecalls to fast5 files
function annotate_raw {
    local fast5path=$1
    local fastqpath=$2
    local summaryFile=$3
    tombo preprocess annotate_raw_with_fastqs --fast5-basedir $fast5path \
          --fastq-filenames <(zcat $fastqpath) \
          --sequencing-summary-filenames $summaryFile \
          --processes 40 \
          --overwrite
}

for i in `ls ../../fast5_input` # i is sampleID: e.g. Leaf_0h_R1
do
    annotate_raw /home/ubuntu/salinity_suppl_analysis/fast5_input/$i /home/ubuntu/salinity_suppl_analysis/fastq_input/$i".fq.gz" /home/ubuntu/salinity_suppl_analysis/fast5_input/$i/sequencing_summary.txt
done

## resquiggle with tombo
function resquiggle {
    local fast5path=$1
    local reference=$2
    tombo resquiggle --overwrite --basecall-group Basecall_1D_001 \
          $fast5path $reference \
          --processes 40 --fit-global-scale --include-event-stdev
}

for i in ../../fast5_input/
do
    resquiggle $i ../stringtie/stringtie.annotated.transcript.fa
done

function nanom6A {
    local fast5path=$1
    local outputpath=$2
    local threads=$3
    find -L $fast5path -name "*.fast5" > $outputpath/files.txt
    extract_raw_and_feature_fast --cpu=$threads \
                                 --fl=$outputpath/files.txt \
                                 -o $outputpath/result \
                                 --clip=10
}

for i in `ls /home/ubuntu/salinity_suppl_analysis/fast5_input/`
do
    if [[ ! -d $i ]]
    then
        mkdir $i
    fi
    nanom6A /home/ubuntu/salinity_suppl_analysis/fast5_input/$i ./$i 40
done

