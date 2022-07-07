#! /usr/bin/env bash

## minimap2 commands for genome and transcriptome mapping

function mappingGenome {
    local sampleName=$(basename $1)
    local outBAM=${sampleName%%.fq.gz}".vsGenome.minimap2.bam"
    minimap2 -ax splice -uf -k14 -t40 --secondary=no -G 20000 --MD \
             /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa \
             $1 |
        samtools view -F 0x4 -u |
        samtools sort -@ 40 -l 9 > $outBAM
}

function mappingTranscriptome {
    local sampleName=$(basename $1)
    local outBAM=${sampleName%%.fq.gz}".vsTranscriptome.minimap2.bam"
    minimap2 -t 40 -ax map-ont -uf --secondary=no \
             /home/ubuntu/data/database/Gmax_508_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.transcript.fa \
             $1 |
        samtools view -F 0x4 -F 0x10 -u |
        samtools sort -@ 40 -l 9 > $outBAM
}

for i in ../../fastq_input/*.fq.gz;
do
    mappingGenome $i
    mappingTranscriptome $i
done
