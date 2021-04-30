#! /bin/bash

## Activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate Nanopore

set -euxo pipefail

POSITIONAL=()

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -m|--module)
        moduleName="$2"
        shift # past argument
        shift # past value
        ;;
    -s|--sample)
        sampleName="$2"
        shift # past argument
        shift # past value
        ;;
    -t|--transcript)
        transFasta="$2"
        shift # past argument
        shift # past value
        ;;
    -r|--ref)
        geneAnno="$2"
        shift # past argument
        shift # past value
        ;;
    -g|--genome)
        genomeFasta="$2"
        shift # past argument
        shift # past value
        ;;
    -n|--threads)
        threads="$2"
        shift # past argument
        shift # past value
        ;;
    -q|--fastq)
        fastqpath="$2"
        shift # past argument
        shift # past value
        ;;
    -o|--output) # Output folder name
        outputpath="$2"
        shift
        shift
        ;;
    -p|--probability) ## probability cutoff for m6A sites
        prob="$2"
        shift
        shift
        ;;
    --support) ## one m6A site supported read number
        supportedReads="$2"
        shift
        shift
        ;;
    *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

fast5path=$1

## moduleName=$1
## sampleName=$2
## fast5path=$3
## fastqpath=$4
## reference_transcript=$5
## threads=$6
## geneAnno=$7
## genomeFasta=$8

if [[ $moduleName == "tombo_preprocess" ]]; then
    tombo preprocess annotate_raw_with_fastqs --fast5-basedir $fast5path \
          --fastq-filenames $fastqpath \
          --overwrite \
          --processes $threads
fi

if [[ $moduleName == "tombo_resquiggle" ]]; then
    tombo resquiggle --overwrite --basecall-group Basecall_1D_000 \
          $fast5path $transFasta \
          --process $threads \
          --fit-global-scale --include-event-stdev
fi

if [[ $moduleName == "nanom6A_extract" ]]; then
    find $fast5path -name "*.fast5" > $outputpath/files.txt
    extract_raw_and_feature_fast --cpu=$threads \
                                 --fl=$outputpath/files.txt \
                                 -o $outputpath/result \
                                 --clip=10
fi

if [[ $moduleName == "nanom6A_predict" ]]; then
    predict_sites --cpu $threads -i $outputpath/result \
                  -o $outputpath/result_final \
                  -r $geneAnno \
                  -g $genomeFasta \
                  --proba $prob \
                  --support $supportedReads
fi
