#! /usr/bin/env Rscript

library(tidyverse)

d1 <- read_tsv("nanoPAS_vs_Gmax_a2v1_three_end.bed",
              col_names = c('chr', 'start', 'end',
                            'name', 'score', 'strand', 'distance')) %>%
    mutate(comparison = "PAS_vs_Annotation")

d2 <- read_tsv("nanoPAS_vs_StringTie_count5_three_end.bed",
               col_names = c('chr', 'start', 'end',
                             'name', 'score', 'strand', 'distance')) %>%
    mutate(comparison = "PAS_vs_StringTie")
d <- rbind(d1,d2) %>%
    select(distance, comparison)

p <- ggplot(d, aes(x=distance, group = comparison,
                   fill = comparison))
p <- p + geom_histogram(bins = 30, alpha =0.6,
                        color = "black",
                        position = "identity")
p <- p + xlim(-500, 500)
p <- p + xlab("Distance (nt)") + ylab("PAS count")
p <- p + scale_y_continuous(labels=scales::comma)
p <- p + scale_fill_brewer(palette = "Set1")
p <- p + theme_bw(base_size=15)
p <- p + theme(panel.grid = element_blank(),
               legend.position = c(0.75,0.85),
               legend.text = element_text(size = 6),
               legend.title = element_blank(),
               axis.text = element_text(color='black'))

ggsave(filename = "nanoPAS_vs_Annotation_three_end.pdf",
       width = 3.5,
       height = 3.5)
