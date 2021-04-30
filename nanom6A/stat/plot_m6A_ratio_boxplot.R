#! /usr/bin/env Rscript

library(tidyverse)

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

## r1 <- ratio_tbl %>%
##     filter(tissue == "C08leaf") %>%
##     select(gene, chr, pos, tissue, sample, ratio) %>%
##     pivot_wider(id_cols = c(gene, chr, pos, tissue), names_from = sample,
##                 values_from = ratio) %>% na.omit() %>%
##     pivot_longer(cols= starts_with("C08leaf"), names_to = "sample",
##                  values_to = "ratio")
## 
## r2 <- ratio_tbl %>%
##     filter(tissue == "C08root") %>%
##     select(gene, chr, pos, tissue, sample, ratio) %>%
##     pivot_wider(id_cols = c(gene, chr, pos, tissue), names_from = sample,
##                 values_from = ratio) %>% na.omit() %>%
##     pivot_longer(cols= starts_with("C08root"), names_to = "sample",
##                  values_to = "ratio")
## 
## ratio_tbl <- rbind(r1,r2)

x <- ratio_tbl %>%
    filter(sample == "C08leaf_0h") %>% pull(ratio)
y <- ratio_tbl %>%
    filter(sample == "C08leaf_1h") %>% pull(ratio)
C08leaf_pvalue <- wilcox.test(x, y)$p.value %>%
                                  signif(digits=3)

x <- ratio_tbl %>%
    filter(sample == "C08root_0h") %>% pull(ratio)
y <- ratio_tbl %>%
    filter(sample == "C08root_1h") %>% pull(ratio)
C08root_pvalue <- wilcox.test(x, y)$p.value %>%
                                  signif(digits=3)

mean_C08leaf_0h <- ratio_tbl %>%
    filter(sample == "C08leaf_0h") %>%
    pull(ratio) %>% mean() %>%
    signif(digits=3)
mean_C08leaf_1h <- ratio_tbl %>%
    filter(sample == "C08leaf_1h") %>%
    pull(ratio) %>% mean() %>%
    signif(digits=3)
mean_C08root_0h <- ratio_tbl %>%
    filter(sample == "C08root_0h") %>%
    pull(ratio) %>% mean() %>%
    signif(digits=3)
mean_C08root_1h <- ratio_tbl %>%
    filter(sample == "C08root_1h") %>%
    pull(ratio) %>% mean() %>%
    signif(digits=3)


p_tbl <- tibble(tissue = c("C08leaf", "C08root"),
                pvalue = c(C08leaf_pvalue, C08root_pvalue),
                mean_0h = c(mean_C08leaf_0h, mean_C08root_0h),
                mean_1h = c(mean_C08leaf_1h, mean_C08root_1h))

count_tbl <- ratio_tbl %>%
    group_by(sample, tissue) %>%
    summarise(count=n())


p <- ggplot(ratio_tbl,
            aes(x=sample, y=ratio,
                fill=sample))
p <- p + geom_violin(alpha=0.8, trim = TRUE)
p <- p + geom_boxplot(alpha=0.8,
                      outlier.shape = NA,
                      width=0.2)
p <- p + xlab("") + ylab("m6A sites ratio")
p <- p + geom_text(data = count_tbl,
                   aes(x=sample, y=1.1,
                       label = paste0("n = ", count)),
                   size = 3)

p <- p + geom_text(data = p_tbl,
                   aes(x=1.5, y=1, fill = NULL,
                       label = paste0("P < ", pvalue)),
                   size = 3)

p <- p + geom_segment(aes(x = 1.15, y = 0.95,
                          xend = 1.85, yend = 0.95))
p <- p + facet_grid(~tissue, scales = "free")

p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               strip.background = element_blank(),
               axis.text = element_text(color = "black"))

ggsave(filename = "m6A_sites_ratio.boxplot.pdf",
       width = 5, height = 3)

p2 <- ggplot(ratio_tbl,
             aes(x=ratio,
                 fill=sample))
p2 <- p2 + geom_density(alpha=0.8)
p2 <- p2 + geom_text(data = p_tbl,
                     aes(x=after_stat(0.8), y=after_stat(1.3), fill = NULL,
                         label = paste0("P < ", pvalue, "\n",
                                        "0h_mean: ", mean_0h, "\n",
                                        "1h_mean: ", mean_1h, "\n")),
                     size = 3)
p2 <- p2 + xlab("Methylation level of m6A sites") + ylab("Density")
p2 <- p2 + facet_grid(~tissue, scales = "free")
p2 <- p2 + theme_bw(base_size=15)
p2 <- p2 + theme(legend.position = "bottom",
                 legend.title = element_blank(),
                 panel.grid = element_blank(),
                 strip.background = element_blank(),
                 axis.text = element_text(color = "black"))

##p2 <- p2 + scale_fill_brewer(palette="Set1")


ggsave(filename = "m6A_sites_ratio.density.pdf",
       width = 6, height = 3.5)
