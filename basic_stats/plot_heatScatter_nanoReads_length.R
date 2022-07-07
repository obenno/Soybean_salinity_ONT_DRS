#! /usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

d <- read_tsv(args[1],
              col_names = c("readID", "transcript", "readLength", "transLength"))
message(nrow(d), " reads included.")

p <- ggplot(d, aes(x=transLength, y=readLength))
##p <- p + stat_density_2d_filled()
p <- p + stat_density_2d(geom = "raster",
                         aes(fill = after_stat(density)),
                         contour = FALSE)
##p <- p + scale_fill_viridis_c()

## Use geom_tile instead
##p <- p + geom_tile(aes(readLength, identity, fill = after_stat(density)))
##p <- p + scale_fill_distiller(palette="Purples", direction = 1)
p <- p + scale_fill_gradient(low="white", high = "#000066")
##p <- p + xlim(c(0,4000)) + ylim(c(0,4000))
p <- p + scale_x_continuous(limits=c(0,4000),
                            labels=scales::comma)
p <- p + scale_y_continuous(limits=c(0,4000),
                            labels=scales::comma)

p <- p + xlab("Expected read length (bases)") + ylab("Observed read length (bases)")

p <- p+theme_bw(base_size=15)+theme(legend.position = "none",
                        panel.grid = element_blank(),
                        panel.border = element_rect(colour = "black", fill=NA, size = 0.6),
                        axis.ticks = element_line(size = 0.6),
                        axis.text = element_text(color = "black"))

ggsave(p, filename = args[2], width = 3.5, height = 3.5)
