#! /usr/bin/env Rscript

library(tidyverse)

psi <- read_tsv("all_events_psi.tsv")

psi_long <- psi %>% pivot_longer(-c(eventID, geneID, AStype, inclusionExon, exclusionExon), names_to = "sample", values_to = "psi")


aggregate_fdr <- function(comparison, samplePattern){
    fdr_tbl <- read_tsv(paste0("rmats_output/", "rMATS_Result_FDR.",
                               comparison, ".txt"))
    outTbl <- psi_long %>% filter(str_detect(sample, samplePattern)) %>%
        right_join(fdr_tbl, by=c("eventID" = "Exon")) %>%
        dplyr::select(-c(IJC1, SJC1, IJC2, SJC2, IncFormLen, SkipFormLen)) %>%
        pivot_wider(c(eventID, geneID, AStype, inclusionExon, exclusionExon, PValue, FDR), names_from = "sample", values_from = "psi")
    return(outTbl)
}

outTbl.Leaf_1h_vs_0h <- aggregate_fdr("Leaf_1h_vs_Leaf_0h", "Leaf_")

outTbl.Root_1h_vs_0h <- aggregate_fdr("Root_1h_vs_Root_0h", "Root_")

outTbl.Root_0h_vs_Leaf_0h <- aggregate_fdr("Root_0h_vs_Leaf_0h", "_0h")

outTbl.Root_1h_vs_Leaf_1h <- aggregate_fdr("Root_1h_vs_Leaf_1h", "_1h")

write_tsv(outTbl.Leaf_1h_vs_0h, "all_events.Leaf_1h_vs_Leaf_0h.FDR.tsv")
write_tsv(outTbl.Root_1h_vs_0h, "all_events.Root_1h_vs_Root_0h.FDR.tsv")
write_tsv(outTbl.Root_0h_vs_Leaf_0h, "all_events.Root_0h_vs_Leaf_0h.FDR.tsv")
write_tsv(outTbl.Root_1h_vs_Leaf_1h, "all_events.Root_1h_vs_Leaf_1h.FDR.tsv")
