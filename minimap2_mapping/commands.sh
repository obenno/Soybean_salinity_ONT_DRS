#! /usr/bin/env bash

## minimap2 commands for genome and transcriptome mapping

minimap2 -ax splice -uf -k14 -t20 --MD /home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa ../fastq/C08leaf_1h.fq.gz | samtools sort -@ 20 -l 9 | samtools view -b > C08leaf_1h.minimap2.bam

minimap2 -ax splice -uf -k14 -t20 --MD /home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa ../fastq/C08leaf_0h.fq.gz | samtools sort -@ 20 -l 9 | samtools view -b > C08leaf_0h.minimap2.bam

minimap2 -ax splice -uf -k14 -t20 --MD /home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa ../fastq/C08root_0h.fq.gz | samtools sort -@ 20 -l 9 | samtools view -b > C08root_0h.minimap2.bam

minimap2 -ax splice -uf -k14 -t20 --MD /home/ubuntu/data/database/Gmax_275_Wm82.a2.v1/assembly/Gmax_275_v2.0.fa ../fastq/C08root_1h.fq.gz | samtools sort -@ 20 -l 9 | samtools view -b > C08root_1h.minimap2.bam

minimap2 -t 30 -ax map-ont ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.transcript.fa ../fastq/C08leaf_0h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_Gmax_a2v1_Transcripts.C08leaf_0h.minimap2.bam

minimap2 -t 30 -ax map-ont ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.transcript.fa ../fastq/C08leaf_1h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_Gmax_a2v1_Transcripts.C08leaf_1h.minimap2.bam

minimap2 -t 30 -ax map-ont ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.transcript.fa ../fastq/C08root_0h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_Gmax_a2v1_Transcripts.C08root_0h.minimap2.bam

minimap2 -t 30 -ax map-ont ~/data/database/Gmax_275_Wm82.a2.v1/annotation/Gmax_275_Wm82.a2.v1.transcript.fa ../fastq/C08root_1h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_Gmax_a2v1_Transcripts.C08root_1h.minimap2.bam

minimap2 -t 30 -ax map-ont stringtie.annotated.transcript.fa ../fastq/C08leaf_0h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_StringTieTranscripts.C08leaf_0h.minimap2.bam

minimap2 -t 30 -ax map-ont stringtie.annotated.transcript.fa ../fastq/C08leaf_1h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_StringTieTranscripts.C08leaf_1h.minimap2.bam

minimap2 -t 30 -ax map-ont stringtie.annotated.transcript.fa ../fastq/C08root_0h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_StringTieTranscripts.C08root_0h.minimap2.bam

minimap2 -t 30 -ax map-ont stringtie.annotated.transcript.fa ../fastq/C08root_1h.fq.gz | samtools view -F 0x4 -F 0x10 -b > NanoReads_vs_StringTieTranscripts.C08root_1h.minimap2.bam
