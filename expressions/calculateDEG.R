#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)
library(DESeq2)

## read tx2gene table, the first two column will be
## used, col names don't matter
tx2gene <- read_tsv("tx2gene.tsv")

## read sample infor
sampleInfo <- data.frame(
    id= list.files(pattern = "^(Leaf|Root)_[01]h_R[123]$")
) %>%
    mutate(
        condition = rep(c("Leaf_0h", "Leaf_1h", "Root_0h", "Root_1h"),
                        each=3)
    )

## salmon output path and sample name
files <- file.path(sampleInfo$id, "quant.sf")
names(files) <- sampleInfo$id

## Generate transcript TPM table
tx.salmon <- tximport(files, type="salmon", txOut = TRUE)
nano_tx_tpm <- tx.salmon$abundance %>%
    as.data.frame() %>%
    rownames_to_column("transcript") %>%
    pivot_longer(cols = -transcript,
                 names_to = "sample", values_to = "nano_TPM")

## Filter nanopore transcript by long read counts >=5
nano_counts5 <- read_tsv(file = "../stringtie/stringtie.transcript.counts_more_than_five.lst", col_names = c("transcript", "count"))

nano_tx_tpm <- nano_tx_tpm %>%
    filter(transcript %in% nano_counts5$transcript)

nano_tx_max_tpm  <- nano_tx_tpm %>%
    group_by(transcript) %>%
    summarize(max_TPM=max(nano_TPM))

nano_tx_maxTPM1 <- nano_tx_max_tpm %>%
    filter(max_TPM>=1) %>%
    select(transcript) %>%
    pull()

nano_tpm_tbl <- nano_tx_tpm %>%
    pivot_wider(names_from = sample,
                values_from = nano_TPM)

nano_tpm_tbl[is.na(nano_tpm_tbl)] <- 0

nano_tpm_tbl_long <- nano_tpm_tbl %>%
    pivot_longer(cols = -transcript,
                 names_to = "sample", values_to = "nano_TPM")

write.table(nano_tx_maxTPM1, file = "Nano_maxTPM1.lst",
            quote = F, sep = "\t",
            row.names = F, col.names = F)

write.table(nano_tpm_tbl, file = "Nano_TPM.tsv",
            quote = F, sep ="\t",
            row.names = F,
            col.names = T)
write.table(nano_tpm_tbl_long,  file = "Nano_TPM.long.tsv",
            quote = F, sep ="\t",
            row.names = F,
            col.names = T)
write.table(tx.salmon$counts, file = "Nano_salmonCounts.tsv",
            quote = F, sep ="\t",
            row.names = T, col.names = T)


## read gene location file
geneLOC <- read_tsv("geneLOC.tsv", col_names = c("geneID", "location"))

## Get gene counts and TPM table
gene.salmon <- tximport(files, type="salmon", tx2gene = tx2gene)
gene_counts <- gene.salmon$counts
gene_tpm <- gene.salmon$abundance

gene2ref <- tx2gene %>%
    dplyr::select(GENEID, Ref_transcript) %>%
    mutate(Ref_transcript = str_replace(Ref_transcript,
                                        regex("\\.[0-9]+\\.Wm82\\.a4\\.v1"),
                                        "\\.Wm82\\.a4\\.v1")) %>%
    distinct() %>% group_by(GENEID) %>%
    summarise(ref = paste(unique(na.omit(Ref_transcript)), collapse = ",")) %>%
    mutate(ref=case_when(str_length(ref)==0 ~ "None", TRUE ~ ref))

gene_tpm <- gene_tpm %>%
    as.data.frame %>%
    rownames_to_column("geneID") %>%
    as_tibble() %>%
    left_join(gene2ref, by=c("geneID" = "GENEID")) %>%
    dplyr::relocate(ref, .after = geneID)

gene_counts <- gene_counts %>%
    as.data.frame %>%
    rownames_to_column("geneID") %>%
    as_tibble() %>%
    left_join(gene2ref, by=c("geneID" = "GENEID")) %>%
    dplyr::relocate(ref, .after = geneID)

write_tsv(gene_tpm, file = "gene.tpm.tsv")

write_tsv(gene_counts, file = "gene.counts.tsv")

## Get transcript counts and TPM table

## Prepare sample meta data for DEseq2
## coldata <- sampleInfo$sample_alias %>%
##     str_split("_", simplify = T) %>%
##     as.data.frame()
## colnames(coldata) <- c("tissue","timepoint", "replicate")
## rownames(coldata) <- sampleInfo$sample_alias
## coldata <- coldata %>%
##     mutate(tissue=str_replace(tissue, "C08", "")) %>%
##     mutate(condition=str_replace(rownames(coldata), "_r[123]", ""))

## Read data into dds object
## use simplese design, only condition
## Test differentially expressed genes (DEGs)
dds <- DESeqDataSetFromTximport(gene.salmon, sampleInfo, ~condition)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

dds$condition <- factor(dds$condition,
                        levels = c("Leaf_0h", "Leaf_1h",
                                   "Root_0h", "Root_1h"))


## Perform DESeq2 DE analysis
dds <- DESeq(dds)

## write output files
deseq_compare <- function(dds, condition="condition",
                          treatment, control, FDR=0.05, LFC=1,
                          gene_tpm = gene_tpm, gene_loc = geneLOC,
                          result_out, up_out, down_out){
    res <- lfcShrink(dds, contrast=c(condition, treatment, control))

    tpm_data <- gene_tpm %>%
        dplyr::select(`geneID`, `ref`,
                      starts_with(c(control, treatment)))

    res <- res %>%
        as.data.frame %>%
        rownames_to_column("geneID") %>%
        as_tibble() %>%
        left_join(tpm_data, by="geneID") %>%
        dplyr::relocate(`ref`, .after = `geneID`) %>%
        left_join(gene_loc, by="geneID") %>%
        dplyr::relocate(`location`, .after = `geneID`)

    write.table(res, file = result_out, quote = F, sep="\t", row.names = F)
    de_up <- res %>%
        dplyr::select(`geneID`, `location`, `ref`,
                      `log2FoldChange`, `pvalue`, `padj`,
                      starts_with(c(control, treatment))) %>%
        dplyr::filter(!is.na(`padj`),
                      `padj` <= FDR,
                      `log2FoldChange` >= LFC)
    de_down <- res %>%
        dplyr::select(`geneID`, `location`, `ref`,
                      `log2FoldChange`, `pvalue`, `padj`,
                      starts_with(c(control, treatment))) %>%
        dplyr::filter(!is.na(`padj`),
                      `padj` <=FDR,
                      `log2FoldChange` <=-LFC)
    write.table(de_up, file = up_out, quote = F, sep="\t", row.names = F)
    write.table(de_down, file = down_out, quote = F, sep="\t", row.names = F)
}

deseq_compare(dds, treatment="Leaf_1h", control="Leaf_0h",
              gene_tpm = gene_tpm, gene_loc = geneLOC,
              result_out="Leaf_1hvs0h.DEGresult.out",
              up_out="Leaf_1hvs0h.gene_up.out",
              down_out="Leaf_1hvs0h.gene_down.out")

deseq_compare(dds, treatment="Root_1h", control="Root_0h",
              gene_tpm = gene_tpm, gene_loc = geneLOC,
              result_out="Root_1hvs0h.DEGresult.out",
              up_out="Root_1hvs0h.gene_up.out",
              down_out="Root_1hvs0h.gene_down.out")

## calculate DTE result
deseq_compare_tx <- function(dds_tx, condition="condition",
                             treatment, control, FDR=0.05, LFC=1,
                             trans_tpm = transTPM, gene_loc = geneLOC,
                             result_out, up_out, down_out){
    res <- lfcShrink(dds_tx, contrast=c(condition, treatment, control))

    tpm_data <- transTPM %>%
        dplyr::select(`transID`, `geneID`, `ref`,
                      starts_with(c(control, treatment)))

    res <- res %>%
        as.data.frame %>%
        rownames_to_column("transID") %>%
        as_tibble() %>%
        left_join(tpm_data, by="transID") %>%
        dplyr::relocate(c(`geneID`,`ref`), .after = `transID`) %>%
        left_join(gene_loc, by="geneID") %>%
        dplyr::relocate(`location`, .after = `geneID`)

    write.table(res, file = result_out, quote = F, sep="\t", row.names = F)
    de_up <- res %>%
        dplyr::select(`transID`, `geneID`, `location`,`ref`,
                      `log2FoldChange`, `pvalue`, `padj`,
                      starts_with(c(control, treatment))) %>%
        dplyr::filter(!is.na(`padj`),
                      `padj` <= FDR,
                      `log2FoldChange` >= LFC)
    de_down <- res %>%
        dplyr::select(`transID`,`geneID`, `location`, `ref`,
                      `log2FoldChange`, `pvalue`, `padj`,
                      starts_with(c(control, treatment))) %>%
        dplyr::filter(!is.na(`padj`),
                      `padj` <=FDR,
                      `log2FoldChange` <=-LFC)
    write.table(de_up, file = up_out, quote = F, sep="\t", row.names = F)
    write.table(de_down, file = down_out, quote = F, sep="\t", row.names = F)
}


tx.salmon <- tximport(files, type="salmon",
                      tx2gene = tx2gene, txOut = TRUE)

transTPM <- tx.salmon$abundance %>%
    as.data.frame() %>%
    rownames_to_column("transID") %>%
    as_tibble() %>%
    left_join(tx2gene, by = c("transID" = "TXNAME")) %>%
    dplyr::relocate(c("GENEID", "Ref_transcript"), .after = `transID`) %>%
    dplyr::rename("geneID" = "GENEID",
                  "ref" = "Ref_transcript")

dds_tx <- DESeqDataSetFromTximport(tx.salmon, sampleInfo, ~condition)
keep <- rowSums(counts(dds_tx)) >= 5
dds_tx <- dds_tx[keep,]

dds_tx$condition <- factor(dds_tx$condition,
                           levels = c("Leaf_0h", "Leaf_1h",
                                      "Root_0h", "Root_1h"))

## Perform DESeq2 DE analysis
dds_tx <- DESeq(dds_tx)

deseq_compare_tx(dds_tx, treatment="Leaf_1h", control="Leaf_0h",
                 trans_tpm = transTPM, gene_loc = geneLOC,
                 result_out="Leaf_1hvs0h.DETresult.out",
                 up_out="Leaf_1hvs0h.trans_up.out",
                 down_out="Leaf_1hvs0h.trans_down.out")

deseq_compare_tx(dds_tx, treatment="Root_1h", control="Root_0h",
                 trans_tpm = transTPM, gene_loc = geneLOC,
                 result_out="Root_1hvs0h.DETresult.out",
                 up_out="Root_1hvs0h.trans_up.out",
                 down_out="Root_1hvs0h.trans_down.out")

