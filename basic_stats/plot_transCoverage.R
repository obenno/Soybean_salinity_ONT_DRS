#! /usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

d <- read_tsv(file = args[1], col_names = c("transcript", "len", "coverage", "source", "lengthGroups"))

d <- d %>%
    mutate(lengthGroups=str_replace(lengthGroups, ",", ", ")) %>%
    mutate(lengthGroups=factor(lengthGroups,
                               levels = c("(0, 500]","(500, 1000]",
                                          "(1000, 1500]","(1500, 2000]",
                                          "(2000, max]")))

## Modify content
d <- d %>%
    mutate(source = str_replace(source,
                                "Annotated", "StringTie (annotated)"),
           source = str_replace(source,
                                "Novel", "StringTie (novel)")) %>%
    mutate(source = factor(source,
                           levels = c("Reference",
                                      "StringTie (annotated)",
                                      "StringTie (novel)")))


count_tbl <- d %>%
    group_by(lengthGroups, source) %>%
    summarise(count=n())

p <- ggplot(d, aes(x=lengthGroups, y=coverage, fill=source))
p <- p + geom_violin(alpha = 0.8)
p <- p + geom_text(data=count_tbl,
                   aes(x = lengthGroups, y = 1.05,
                       label=scales::comma(count, accuracy=1)),
                   color="black")

p <- p + facet_grid(cols=vars(source))

p <- p + xlab("Transcript length")
p <- p + ylab("Transcript coverage per alignment")
p <- p + theme_bw(base_size=15)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               panel.border = element_rect(colour = "black", fill=NA, size = 0.6),
               axis.ticks = element_line(size = 0.6),
               axis.text = element_text(color = "black"),
               axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
p <- p + theme(strip.background = element_blank())
p <- p + scale_fill_brewer(palette = "Dark2")
##p <- p + scale_fill_brewer(palette = "Dark2")

ggsave(filename = args[2], width = 9, height = 4.5)

