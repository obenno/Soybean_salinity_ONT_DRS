#! /usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

d <- read_tsv(args[1], col_names = F)
colnames(d) <- c("readID", "identity", "readLength")

p <- ggplot(d, aes(x=readLength, y=identity))
##p <- p + stat_density_2d_filled()
p <- p + stat_density_2d(geom = "raster",
                         aes(fill = after_stat(density)),
                         contour = FALSE)
##p <- p + scale_fill_viridis_c()

## Use geom_tile instead
##p <- p + geom_tile(aes(readLength, identity, fill = after_stat(density)))
##p <- p + scale_fill_distiller(palette="Purples", direction = 1)
p <- p + scale_fill_gradient(low="white", high = "#000066")
p <- p + ylim(c(0.5,1))
p <- p + scale_x_continuous(limits=c(0,3000),
                            ##labels=scales::label_number(scale=1/1000, suffix="K"))
                            labels=scales::comma)
p <- p + xlab("Read length (bases)") + ylab("RNA read identity")

p <- p+theme_bw(base_size=15)+theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.border = element_rect(colour = "black", fill=NA, size = 0.6),
                        axis.ticks = element_line(size = 0.6),
                        axis.text = element_text(color = "black"))

ggsave(p, filename = args[2], width = 3.5, height = 3.5)
