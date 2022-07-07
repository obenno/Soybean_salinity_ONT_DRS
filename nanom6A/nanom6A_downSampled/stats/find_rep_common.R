#! /usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

rep1_name <- args[1]
rep2_name <- args[2]
rep3_name <- args[3]

rep1_ratioFile <- paste0("m6A_ratio.", rep1_name, ".tsv")
rep2_ratioFile <- paste0("m6A_ratio.", rep2_name, ".tsv")
rep3_ratioFile <- paste0("m6A_ratio.", rep3_name, ".tsv")

sample1 <- read_tsv(rep1_ratioFile, col_names = c("geneID", "chr", "pos", "motif", "ratio", "mod_n", "total_n"))
sample2 <- read_tsv(rep2_ratioFile, col_names = c("geneID", "chr", "pos", "motif", "ratio", "mod_n", "total_n"))
sample3 <- read_tsv(rep3_ratioFile, col_names = c("geneID", "chr", "pos", "motif", "ratio", "mod_n", "total_n"))

meth_combined <- full_join(
    sample1, sample2,
    by=c("geneID", "chr", "pos", "motif"),
    suffix = c(".r1", ".r2")
) %>%
    full_join(
        sample3,
        by=c("geneID", "chr", "pos", "motif")
    ) %>%
    dplyr::rename(
               "ratio.r3" = "ratio",
               "mod_n.r3" = "mod_n",
               "total_n.r3" = "total_n"
           ) %>% na.omit()

write.table(meth_combined, file = args[4], sep="\t",
            quote = FALSE, row.names = FALSE)
