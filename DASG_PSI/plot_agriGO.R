#! /usr/bin/env Rscript

library(tidyverse)

suppressPackageStartupMessages(library(optparse))

args <- commandArgs(trailingOnly = T)
option.list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input agri GO files, separate by comma"),
    make_option(c("-n", "--names"), type="character", default=NULL,
                help="category names, separatre by comma"),
    make_option(c("-c", "--category"), type="character", default=NULL,
                help="GO category, BP, MF, CC or All"),
    ## deprecate simplify
    ##make_option(c("-s", "--simplify"), action="store_true", default=FALSE,
    ##            help="If perform simplify [%default]"),
    make_option(c("-f","--flip"), action="store_true", default=FALSE,
                help="Whether flip the plot [%default]"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output file name"),
    make_option(c("-w", "--width"), type="double", default=9,
                help="output file width"),
    make_option(c("-e", "--height"), type="double", default=6,
                help="output file height")
)
usage <- "usage: %prog -i A_GO.txt,B_GO.txt -n TypeA,TypeB -o output.pdf\n"
parser <- OptionParser(option_list=option.list,
                       description="Comparative GO enrichment",
                       usage=usage)
opt <- parse_args(parser, args=args, positional_arguments=T)

inputFiles <- unlist(strsplit(opt$options$input, ",",
                       fixed = TRUE))
CatNames <- unlist(strsplit(opt$options$names, ",",
                     fixed = TRUE))

if(length(inputFiles)!=length(CatNames)){
    stop("The category names doesn't match input list length...")
}

if(is.null(opt$options$category)){
    message("GO category not indicated, use all three")
    opt$options$category <- "All"
}else if(!(opt$options$category %in% c("BP", "MF", "CC", "All"))){
    stop("Please indicate GO category")
}

A_GO <- read_tsv(inputFiles[1]) %>%
    mutate(sample=CatNames[1])

B_GO <- read_tsv(inputFiles[2]) %>%
    mutate(sample=CatNames[2])

GO_result <- rbind(A_GO, B_GO) %>%
    filter(FDR<=0.05) %>%
    select(sample, GO_acc, Term, term_type, FDR, queryitem) %>%
    mutate(FDR = -log10(FDR))

if(opt$options$category == "BP"){
    GO_result <- GO_result %>%
        filter(term_type == "P")
}else if(opt$options$category == "MF"){
    GO_result <- GO_result %>%
        filter(term_type == "F")
}else if(opt$options$category == "CC"){
    GO_result <- GO_result %>%
        filter(term_type == "C")
}

GO_result <- GO_result %>%
    arrange(sample, queryitem) %>%
    mutate(Term = factor(Term, levels = unique(Term)))

##p <- ggplot(GO_result, aes(x=Term, y=FDR, fill=sample))
##p <- p + geom_bar(stat = "identity",
##                  position = "dodge")
##p <- p + geom_hline(yintercept = -log10(0.05))
p <- ggplot(GO_result, aes(y = sample,
                           x = Term,
                           size = queryitem,
                           fill = FDR))
p <- p + geom_point(shape = 21)
p <- p + scale_fill_distiller(palette = 'Blues',
                              name = '-log10(FDR)')
p <- p + scale_size_binned(name = 'NumGenes')
##continuous(breaks=c(5, 100, 1000, 2000, 3000))

if(opt$options$category == "All"){
    p <- p + facet_grid(~term_type,
                        scales = "free",
                        space = "free")
}
p <- p + xlab("GO_terms [P]") + ylab("Category")
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank(),
               axis.text.x = element_text(angle=90,
                                          vjust = 0.5,
                                          hjust =1),
               axis.text = element_text(color = 'black'))

p <- p + theme(legend.position = c(0.8, -4),
               legend.box.just = "bottom",
               legend.direction = "horizontal")
ggsave(file = opt$options$output,
       width = opt$options$width,
       height = opt$options$height)
