#! /usr/bin/env Rscript

library(tidyverse)

d <- read_tsv("fusion_candidates.rawCounts.Long.tsv")

d  <- d %>%
    select(transID, refGene, counts, sampleID) %>%
    pivot_wider(id_cols = c(transID, refGene),
                names_from = sampleID,
                values_from = counts)

## replace NA to 0
d[is.na(d)] <- 0

write_tsv(d, file = "fusion_candidates.rawCounts.tsv")

