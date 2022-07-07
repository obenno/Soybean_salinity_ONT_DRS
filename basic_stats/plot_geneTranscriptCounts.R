#! /usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

d <- read_tsv(args[1], col_names = c("geneID", "transCount", "group"))

d <- d %>%
    mutate(group = str_replace(group, "Annotated", "StringTie (annotated)"),
           group = str_replace(group, "Novel", "StringTie (novel)")) %>%
    mutate(group=factor(group,
                        levels = c("Reference",
                                   "StringTie (annotated)",
                                   "StringTie (novel)")))

p <- ggplot(d, aes(x=transCount, group=group, fill=group))
p <- p + geom_bar(position = "dodge2",
                  alpha=0.8, color = "black",
                  size = 0.4)

p <- p + scale_x_continuous(breaks = c(1:10), limits = c(0.5,10))
p <- p + scale_y_continuous(labels = scales::comma)
p <- p + theme_bw(base_size=15)
p <- p + theme(legend.position = c(0.65,0.8),
               legend.title = element_blank(),
               legend.text = element_text(size=8),
               panel.grid = element_blank(),
               panel.border = element_rect(colour = "black",
                                           fill=NA, size = 0.6),
               axis.ticks = element_line(size = 0.6),
               axis.text = element_text(color = "black"))


p <- p + scale_fill_brewer(palette = "Dark2")
p <- p + xlab("Transcript counts for genes")
p <- p + ylab("Gene counts")

ggsave(args[2], width = 3.5, height = 3.5)

