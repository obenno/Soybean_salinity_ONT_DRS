#! /usr/bin/env Rscript

library(tidyverse)

expr_tbl <- read_tsv("Nanopore_salmonminimap2_TPM.long.tsv",
                     col_names = c("gene", "transcript",
                                   "sample", "TPM"))
geneExpr <- expr_tbl %>%
    group_by(gene, sample) %>%
    summarize(TPM=sum(TPM))

readData <- function(fileIn, sample){
    d <- read_tsv(fileIn,
                  col_names = c("gene", "chr", "pos", "motif", "ratio",
                                "methyl_count", "total_count")) %>%
        mutate(sample = sample)
    return(d)
}

sample1 <- readData("m6A_ratio.C08leaf_0h.tsv", "C08leaf_0h")
sample2 <- readData("m6A_ratio.C08leaf_1h.tsv", "C08leaf_1h")
sample3 <- readData("m6A_ratio.C08root_0h.tsv", "C08root_0h")
sample4 <- readData("m6A_ratio.C08root_1h.tsv", "C08root_1h")

ratio_tbl <- rbind(sample1, sample2,
                   sample3, sample4)
ratio_tbl <- ratio_tbl %>%
    mutate(tissue=str_sub(sample, 1, 7))

ratio_geneExpr <- ratio_tbl %>%
    select(gene, sample, ratio) %>%
    group_by(gene, sample) %>%
    summarize(ratio = mean(ratio)) %>%
    left_join(y=geneExpr, by=c("gene", "sample"))

set.seed(123)

cor_tbl <- ratio_geneExpr %>%
    group_by(sample) %>%
    dplyr::summarize(Cor=cor(ratio, TPM,
                             method = "spearman",
                             use = "complete.obs")) %>%
    mutate(Cor=signif(Cor, digits = 2))

p <- ggplot(ratio_geneExpr,
            aes(x=ratio,
                y=log2(TPM+1)))
p <- p + geom_point(alpha=0.5, color = "#4682B4")
p <- p + geom_text(data=cor_tbl,
                   aes(x=0.7, y=15,
                       label=paste0("Spearman cor = ", Cor)))
p <- p + xlab("Methylation level of m6A sites")
p <- p + facet_grid(~sample)

p <- p + theme_bw(base_size = 15)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               panel.border = element_rect(colour = "black",
                                           fill=NA, size = 0.6),
               axis.ticks = element_line(size = 0.6),
               axis.text = element_text(color = "black"),
               strip.background = element_blank())

ggsave(filename = "m6A_ratio_vs_geneExpr.pdf",
       width = 12, height = 3.5)
