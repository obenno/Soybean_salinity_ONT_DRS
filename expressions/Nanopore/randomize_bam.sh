#! /usr/bin/env bash

inputBAM=$1
outBAM=$2

header=$(mktemp -p ./)
alignments=$(mktemp -p ./)

samtools view -H $inputBAM > $header
samtools view $inputBAM |shuf > $alignments
## Add filter here
cat $header $alignments | samtools view -F 0x4 -F 0x10 -F 0x800 -F 0x100 -b > $outBAM
rm $header $alignments
