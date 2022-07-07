#! /usr/bin/env Rscript

library(tidyverse)
library(foreach)

expr_tbl <- read_tsv("/home/ubuntu/salinity_suppl_analysis/analysis/expressions/Nano_TPM.withRefGeneID.tsv") %>%
        pivot_longer(-c(RefGene, GeneID, transcript),
                 names_to = "sample",
                 values_to = "TPM")

geneExpr <- expr_tbl %>%
    group_by(RefGene, GeneID, sample) %>%
    summarize(TPM=sum(TPM))

m6A_outFiles <- dir(pattern = "m6A_ratio.*.tsv")

ratio_tbl <- foreach(i=1:length(m6A_outFiles), .combine="bind_rows") %do% {
    sample <- unlist(str_split(m6A_outFiles[i], fixed(".")))[2]
    read_tsv(m6A_outFiles[i],
             col_names = c("gene", "chr", "pos", "motif", "ratio",
                           "methyl_count", "total_count")) %>%
        mutate(sample = sample)
}

## multiple m6A sites located in the same gene were averaged
ratio_geneExpr <- ratio_tbl %>%
    select(gene, sample, ratio) %>%
    group_by(gene, sample) %>%
    summarize(ratio = mean(ratio)) %>%
    inner_join(y=geneExpr, by=c("gene" = "GeneID", "sample"))

set.seed(123)

cor_tbl <- ratio_geneExpr %>%
    group_by(sample) %>%
    dplyr::summarize(Cor=cor(ratio, log2(TPM+1),
                             method = "spearman",
                             use = "complete.obs")) %>%
    mutate(Cor=signif(Cor, digits = 2))

## Add sampleName and replciate infor for facet plotting
ratio_geneExpr <- ratio_geneExpr %>%
    mutate(sampleName=str_replace(sample, "_R[123]$", ""),
           replicate=str_extract(sample, "R[123]$"))

cor_tbl <- cor_tbl %>%
    mutate(sampleName=str_replace(sample, "_R[123]$", ""),
           replicate=str_extract(sample, "R[123]$"))

p <- ggplot(ratio_geneExpr,
            aes(x=ratio,
                y=log2(TPM+1)))
##p <- p + geom_point(alpha=0.5, color = "#4682B4")
p <- p + geom_point(alpha=0.5, color = "#7570b3")
p <- p + geom_text(data=cor_tbl,
                   aes(x=0.7, y=15,
                       label=paste0("Spearman cor = ", Cor)))
p <- p + xlab("Methylation level of m6A sites")
p <- p + facet_grid(sampleName~replicate)

p <- p + theme_bw(base_size = 15)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               panel.border = element_rect(colour = "black",
                                           fill=NA, size = 0.6),
               axis.ticks = element_line(size = 0.6),
               axis.text = element_text(color = "black"),
               strip.background = element_blank())

ggsave(filename = "m6A_ratio_vs_geneExpr.nanom6A.png",
       width = 9, height = 12)
