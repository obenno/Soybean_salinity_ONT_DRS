#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)
library(DESeq2)

## read tx2gene table, the first two column will be
## used, col names don't matter
tx2gene <- read_tsv("tx2gene.tsv")

## read sample infor
sampleInfo <- read_tsv(file = "filereport_read_run_PRJNA432861_tsv.txt")
sampleInfo <- sampleInfo %>%
    mutate(sample_alias=str_replace(sample_alias,
                                    "_Salt_Transcriptome_", ""),
           sample_alias=str_replace(sample_alias, "Leaf", "leaf"),
           sample_alias=str_replace(sample_alias,"Root", "root"))

## salmon output path and sample name
files <- file.path(sampleInfo$run_accession, "quant.sf")
names(files) <- sampleInfo$sample_alias

## Get gene counts and TPM table
gene.salmon <- tximport(files, type="salmon", tx2gene = tx2gene)
gene_counts <- gene.salmon$counts
gene_tpm <- gene.salmon$abundance
write.table(gene_tpm, file = "gene.tpm.tsv",
            quote =F, row.names = T, col.names = T,
            sep = "\t")
write.table(gene_counts, file = "gene.counts.tsv",
            quote =F, row.names = T, col.names = T,
            sep = "\t")

## Get transcript counts and TPM table
tx.salmon <- tximport(files, type="salmon",
                      tx2gene = tx2gene, txOut = TRUE)
tx_counts <- tx.salmon$counts %>% as.data.frame()
tx_tpm <- tx.salmon$abundance %>% as.data.frame()
write.table(tx_tpm, file = "transcript.tpm.tsv",
            quote =F, row.names = T, col.names = T,
            sep = "\t")
write.table(tx_counts, file = "transcript.counts.tsv",
            quote =F, row.names = T, col.names = T,
            sep = "\t")

## Prepare sample meta data for DEseq2
coldata <- sampleInfo$sample_alias %>%
    str_split("_", simplify = T) %>%
    as.data.frame()
colnames(coldata) <- c("tissue","timepoint", "replicate")
rownames(coldata) <- sampleInfo$sample_alias
coldata <- coldata %>%
    mutate(tissue=str_replace(tissue, "C08", "")) %>%
    mutate(condition=str_replace(rownames(coldata), "_r[123]", ""))

## Read data into dds object
## use simplese design, only condition
## Test differentially expressed genes (DEGs)
dds <- DESeqDataSetFromTximport(gene.salmon, coldata, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

deseq_compare <- function(dds, condition="condition",
                          treatment, control, FDR=0.01, LFC=1,
                          result_out,up_out, down_out){
    res <- results(dds, contrast=c(condition, treatment, control))
    write.table(res, file = result_out, quote = F, sep="\t", row.names = T)
    de_up <- subset(res, padj!= "NA" & padj<=FDR & log2FoldChange>=LFC)
    de_down <- subset(res, padj!= "NA" & padj<=FDR & log2FoldChange<=-LFC)
    write.table(de_up, file = up_out, quote = F, sep="\t", row.names = T)
    write.table(de_down, file = down_out, quote = F, sep="\t", row.names = T)
}

deseq_compare(dds, treatment="C08leaf_1h", control="C08leaf_0h",
              result_out="C08leaf_1hvs0h.DEGresult.out",
              up_out="C08leaf_1hvs0h.gene_up.out",
              down_out="C08leaf_1hvs0h.gene_down.out")
deseq_compare(dds, treatment="C08root_1h", control="C08root_0h",
              result_out="C08root_1hvs0h.DEGresult.out",
              up_out="C08root_1hvs0h.gene_up.out",
              down_out="C08root_1hvs0h.gene_down.out")

## Test differentially expressed isoforms (DEI)
dds <- DESeqDataSetFromTximport(tx.salmon, coldata, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
deseq_compare(dds, treatment="C08leaf_1h", control="C08leaf_0h",
              result_out="C08leaf_1hvs0h.DEIresult.out",
              up_out = "C08leaf_1hvs0h.transcript_up.out",
              down_out = "C08leaf_1hvs0h.transcript_down.out")
deseq_compare(dds, treatment="C08root_1h", control="C08root_0h",
              result_out="C08root_1hvs0h.DEIresult.out",
              up_out = "C08root_1hvs0h.transcript_up.out",
              down_out = "C08root_1hvs0h.transcript_down.out")

