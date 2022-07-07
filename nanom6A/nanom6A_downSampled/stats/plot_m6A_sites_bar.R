#! /usr/bin/env Rscript

library(tidyverse)

Leaf_0h <- read_tsv("Leaf_0h.ratioCommon.tsv") %>%
    rowwise() %>%
    dplyr::filter(min(total_n.r1,total_n.r2,total_n.r3)>=20) %>%
    select(geneID, chr, pos, motif) %>%
    mutate(sample="Leaf_0h")

Leaf_1h <- read_tsv("Leaf_1h.ratioCommon.tsv") %>%
    rowwise() %>%
    dplyr::filter(min(total_n.r1,total_n.r2,total_n.r3)>=20) %>%
    select(geneID, chr, pos, motif) %>%
    mutate(sample="Leaf_1h")

Root_0h <- read_tsv("Root_0h.ratioCommon.tsv") %>%
    rowwise() %>%
    filter(min(total_n.r1,total_n.r2,total_n.r3)>=20) %>%
    select(geneID, chr, pos, motif) %>%
    mutate(sample="Root_0h")

Root_1h <- read_tsv("Root_1h.ratioCommon.tsv") %>%
    rowwise() %>%
    filter(min(total_n.r1,total_n.r2,total_n.r3)>=20) %>%
    select(geneID, chr, pos, motif) %>%
    mutate(sample="Root_1h")


d  <- bind_rows(Leaf_0h, Leaf_1h, Root_0h, Root_1h)

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
               legend.position = "bottom")


ggsave(filename = "m6A_sites_count.bar.pdf",
       width = 6, height = 6)

Leaf_0h_sites <- d %>%
    filter(sample == "Leaf_0h") %>%
    mutate(siteID=paste(geneID, chr, pos, motif, sep="::")) %>%
    pull(siteID)

write.table(Leaf_0h_sites, file = "Leaf_0h_sites.lst",
            sep="\t", row.names = F, col.names = F, quote = F)

Leaf_1h_sites <- d %>%
    filter(sample == "Leaf_1h") %>%
    mutate(siteID=paste(geneID, chr, pos, motif, sep="::")) %>%
    pull(siteID)

write.table(Leaf_1h_sites, file = "Leaf_1h_sites.lst",
            sep="\t", row.names = F, col.names = F, quote = F)

Root_0h_sites <- d %>%
    filter(sample == "Root_0h") %>%
    mutate(siteID=paste(geneID, chr, pos, motif, sep="::")) %>%
    pull(siteID)

write.table(Root_0h_sites, file = "Root_0h_sites.lst",
            sep="\t", row.names = F, col.names = F, quote = F)

Root_1h_sites <- d %>%
    filter(sample == "Root_1h") %>%
    mutate(siteID=paste(geneID, chr, pos, motif, sep="::")) %>%
    pull(siteID)

write.table(Root_1h_sites, file = "Root_1h_sites.lst",
            sep="\t", row.names = F, col.names = F, quote = F)
