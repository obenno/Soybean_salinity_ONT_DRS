#! /usr/bin/env Rscript

library(tidyverse)
library(DRIMSeq)

## This script is to calculate PSI and format count table for rmats
## since PAS framework was replaced by QuantifyPolyA, PAS counts
## will not be included
expr_long <- read_tsv("../expressions/Nano_TPM.long.tsv")

count_long <- read.delim2("../expressions/Nano_salmonCounts.tsv", header = TRUE, stringsAsFactors = FALSE) %>%
    rownames_to_column("transcript") %>%
    pivot_longer(-transcript, names_to = "sample", values_to = "count") %>%
    mutate(count = as.numeric(count))

AS_info <- read_tsv(
    "stringtie.count_5.maxTPM_1.AS.withChrStrand.out",
    col_names = c("geneID", "chr", "strand",
                  "AStype",
                  "inclusionStart", "inclusionEnd",
                  "exclusionStart", "exclusionEnd",
                  "inclusionTrans", "exclusionTrans")
) %>%
    mutate(
        inclusionExon = paste0(
            chr, ":", inclusionStart, "-", inclusionEnd, ":", strand
        ),
        exclusionExon = paste0(
            chr, ":", exclusionStart, "-", exclusionEnd, ":", strand
        )
    ) %>%
    dplyr::select(-c(chr, strand,
                     inclusionStart, inclusionEnd,
                     exclusionStart, exclusionEnd))


AS_inclusionExpr <- AS_info %>%
    mutate(inclusionTrans=str_split(inclusionTrans, ",")) %>%
    dplyr::select(geneID, AStype, inclusionExon, exclusionExon, inclusionTrans) %>%
    unnest(inclusionTrans) %>%
    left_join(expr_long, by=c("inclusionTrans" = "transcript")) %>%
    group_by(geneID, AStype, inclusionExon, exclusionExon, sample) %>%
    summarise(tpm=sum(nano_TPM, na.rm = TRUE)) %>%
    mutate(exonType = "inclusion")

AS_exclusionExpr <- AS_info %>%
    mutate(exclusionTrans=str_split(exclusionTrans, ",")) %>%
    dplyr::select(geneID, AStype, inclusionExon, exclusionExon, exclusionTrans) %>%
    unnest(exclusionTrans) %>%
    left_join(expr_long, by=c("exclusionTrans" = "transcript")) %>%
    group_by(geneID, AStype, inclusionExon, exclusionExon, sample) %>%
    summarise(tpm=sum(nano_TPM, na.rm = TRUE)) %>%
    mutate(exonType = "exclusion")

AS_event_totalTPM <- bind_rows(AS_inclusionExpr, AS_exclusionExpr) %>%
    ungroup() %>%
    group_by(geneID, AStype, inclusionExon, exclusionExon, sample) %>%
    summarise(totalTPM=sum(tpm, na.rm = TRUE))

AS_event_PSI <- AS_inclusionExpr %>%
    left_join(
        AS_event_totalTPM,
        by=c("geneID", "AStype", "inclusionExon", "exclusionExon", "sample")
    ) %>%
    mutate(psi=tpm/totalTPM) %>%
    ungroup() %>%
    pivot_wider(c(geneID, AStype, inclusionExon, exclusionExon),
                names_from = "sample",
                values_from = "psi") %>%
    mutate(eventID = paste0("ASevent_", row_number())) %>%
    dplyr::relocate("eventID", .before = 1)

AS_inclusionCount <- AS_info %>%
    mutate(inclusionTrans=str_split(inclusionTrans, ",")) %>%
    dplyr::select(geneID, AStype,
           inclusionExon, exclusionExon, inclusionTrans) %>%
    unnest(inclusionTrans) %>%
    left_join(count_long, by=c("inclusionTrans" = "transcript")) %>%
    group_by(geneID, AStype, inclusionExon, exclusionExon, sample) %>%
    summarise(count=sum(count, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(id=paste(geneID, AStype, inclusionExon, exclusionExon,
                    sep="|")) %>%
    dplyr::select(id, sample, count) %>%
    pivot_wider(id, names_from = "sample", values_from = "count")

AS_exclusionCount <- AS_info %>%
    mutate(exclusionTrans=str_split(exclusionTrans, ",")) %>%
    dplyr::select(geneID, AStype,
           inclusionExon, exclusionExon, exclusionTrans) %>%
    unnest(exclusionTrans) %>%
    left_join(count_long, by=c("exclusionTrans" = "transcript")) %>%
    group_by(geneID, AStype, inclusionExon, exclusionExon, sample) %>%
    summarise(count=sum(count, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(id=paste(geneID, AStype, inclusionExon, exclusionExon,
                    sep="|")) %>%
    dplyr::select(id, sample, count) %>%
    pivot_wider(id, names_from = "sample", values_from = "count")


inclusionInputData <- AS_inclusionCount %>% arrange(id) %>%
    mutate(
        eventID=paste0("ASevent_", 1:nrow(AS_inclusionCount)),
        exonID=paste0(eventID, "|", "inclusion")
    ) %>%
    dplyr::relocate(c(eventID, exonID)) %>%
    dplyr::select(-id)

exclusionInputData <- AS_exclusionCount %>% arrange(id) %>%
    mutate(
        eventID=paste0("ASevent_", 1:nrow(AS_exclusionCount)),
        exonID=paste0(eventID, "|", "exclusion")
    ) %>%
    dplyr::relocate(c(eventID, exonID)) %>%
    dplyr::select(-id)


DRIMcounts <- bind_rows(
    inclusionInputData, exclusionInputData##,
    ##PAS_proximal_inputData, PAS_distal_inputData
) %>%
    dplyr::rename("gene_id" = `eventID`, "feature_id" = `exonID`) %>%
    as.data.frame

## Force na value to 0
DRIMcounts <- DRIMcounts %>%
    mutate(
        across(
            .cols = everything(),
            .fns = function(x){ifelse(is.na(x), 0 , x)}
        )
    )

set.seed(1234)

sampleNames <- c("Leaf_0h", "Leaf_1h", "Root_0h", "Root_1h")

replicates <- c("R1", "R2", "R3")

sampleInfo <- data.frame(sample_id=paste(sampleNames,
                                  rep(replicates, each=4),
                                  sep="_")) %>%
    mutate(condition = str_replace(sample_id, "_R[123]", "")) %>%
    arrange(condition)

d <- dmDSdata(counts=DRIMcounts, samples=sampleInfo)

n <- 12
n.small <- 3
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=3,
              min_samps_feature_prop=n.small, min_feature_prop=0.05,
              min_samps_gene_expr=n.small, min_gene_expr=5)


countTbl <- counts(d) %>%
    mutate(exonType = str_replace(feature_id, regex("^(ASevent_|aPAS_)+\\d+\\|"), "")) %>%
    dplyr::select(-feature_id) %>%
    pivot_longer(-c(gene_id, exonType),
                 names_to = "library",
                 values_to = "count") %>%
    mutate(sample = str_replace(library, "_R[123]", ""))

rmats_input.Leaf_1h_vs_Leaf_0h <- countTbl %>%
    dplyr::filter(sample %in% c("Leaf_0h", "Leaf_1h")) %>%
    group_by(gene_id, exonType, sample) %>%
    summarise(count=paste(count, collapse=",")) %>%
    mutate(exonType = str_replace(exonType, "distal", "inclusion"),
           exonType = str_replace(exonType, "proximal", "exclusion")) %>%
    pivot_wider(gene_id,
                names_from = c("exonType", "sample"),
                values_from = count) %>%
    mutate(IncFormLen = 1, SkipFormLen = 1) %>%
    dplyr::rename("Exon" = gene_id,
                  "IJC1" = inclusion_Leaf_0h,
                  "SJC1" = exclusion_Leaf_0h,
                  "IJC2" = inclusion_Leaf_1h,
                  "SJC2" = exclusion_Leaf_1h) %>%
    dplyr::select(Exon, IJC1, SJC1, IJC2, SJC2, IncFormLen, SkipFormLen)

rmats_input.Root_1h_vs_Root_0h <- countTbl %>%
    dplyr::filter(sample %in% c("Root_0h", "Root_1h")) %>%
    group_by(gene_id, exonType, sample) %>%
    summarise(count=paste(count, collapse=",")) %>%
    mutate(exonType = str_replace(exonType, "distal", "inclusion"),
           exonType = str_replace(exonType, "proximal", "exclusion")) %>%
    pivot_wider(gene_id,
                names_from = c("exonType", "sample"),
                values_from = count) %>%
    mutate(IncFormLen = 1, SkipFormLen = 1) %>%
    dplyr::rename("Exon" = gene_id,
                  "IJC1" = inclusion_Root_0h,
                  "SJC1" = exclusion_Root_0h,
                  "IJC2" = inclusion_Root_1h,
                  "SJC2" = exclusion_Root_1h) %>%
    dplyr::select(Exon, IJC1, SJC1, IJC2, SJC2, IncFormLen, SkipFormLen)

rmats_input.Root_0h_vs_Leaf_0h <- countTbl %>%
    dplyr::filter(sample %in% c("Leaf_0h", "Root_0h")) %>%
    group_by(gene_id, exonType, sample) %>%
    summarise(count=paste(count, collapse=",")) %>%
    mutate(exonType = str_replace(exonType, "distal", "inclusion"),
           exonType = str_replace(exonType, "proximal", "exclusion")) %>%
    pivot_wider(gene_id,
                names_from = c("exonType", "sample"),
                values_from = count) %>%
    mutate(IncFormLen = 1, SkipFormLen = 1) %>%
    dplyr::rename("Exon" = gene_id,
                  "IJC1" = inclusion_Leaf_0h,
                  "SJC1" = exclusion_Leaf_0h,
                  "IJC2" = inclusion_Root_0h,
                  "SJC2" = exclusion_Root_0h) %>%
    dplyr::select(Exon, IJC1, SJC1, IJC2, SJC2, IncFormLen, SkipFormLen)

rmats_input.Root_1h_vs_Leaf_1h <- countTbl %>%
    dplyr::filter(sample %in% c("Leaf_1h", "Root_1h")) %>%
    group_by(gene_id, exonType, sample) %>%
    summarise(count=paste(count, collapse=",")) %>%
    mutate(exonType = str_replace(exonType, "distal", "inclusion"),
           exonType = str_replace(exonType, "proximal", "exclusion")) %>%
    pivot_wider(gene_id,
                names_from = c("exonType", "sample"),
                values_from = count) %>%
    mutate(IncFormLen = 1, SkipFormLen = 1) %>%
    dplyr::rename("Exon" = gene_id,
                  "IJC1" = inclusion_Leaf_1h,
                  "SJC1" = exclusion_Leaf_1h,
                  "IJC2" = inclusion_Root_1h,
                  "SJC2" = exclusion_Root_1h) %>%
    dplyr::select(Exon, IJC1, SJC1, IJC2, SJC2, IncFormLen, SkipFormLen)

write_tsv(rmats_input.Leaf_1h_vs_Leaf_0h, "rmats_input.Leaf_1h_vs_Leaf_0h.tsv")

write_tsv(rmats_input.Root_1h_vs_Root_0h, "rmats_input.Root_1h_vs_Root_0h.tsv")

write_tsv(rmats_input.Root_0h_vs_Leaf_0h, "rmats_input.Root_0h_vs_Leaf_0h.tsv")

write_tsv(rmats_input.Root_1h_vs_Leaf_1h, "rmats_input.Root_1h_vs_Leaf_1h.tsv")

write_tsv(AS_event_PSI, "AS_event_PSI.tsv")

