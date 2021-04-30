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

d <- rbind(sample1, sample2,
           sample3, sample4)

p <- ggplot(d,
            aes(x=motif, group=sample,
                fill=sample))

p <- p + geom_bar(position = "dodge",
                  color = "black",
                  alpha = 0.8)
p <- p + scale_fill_brewer(palette="Dark2")

p <- p + theme_bw(base_size=15)
p <- p + xlab("m6A motifs") + ylab("Counts")
p <- p + scale_y_continuous(labels=scales::comma)
p <- p + theme(axis.text = element_text(color = "black"),
               axis.text.x = element_text(angle = 90,
                                           vjust = 0.5,
                                           hjust =1),
               panel.grid = element_blank(),
               legend.title = element_blank(),
               legend.position = c(0.85, 0.85))


ggsave(filename = "m6A_sites_count.bar.pdf",
       width = 5, height = 5)
