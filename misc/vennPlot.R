#! /usr/env/Rscript

## This script is to generate a venn diagram for up to 5 lists
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(optparse))
library(ggplot2)
## suppress log file
## flog.threshold(ERROR)

## make options
args <- commandArgs(trailingOnly = T)
option.list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="input list files, separate by comma"),
    make_option(c("-n", "--names"), type="character", default=NULL,
                help="category names, separatre by comma"),
    make_option(c("-t", "--title"), type="character", default=NULL,
                help="title of the venn diagram"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output file name"),
    make_option(c("--width"), type="numeric", default=1.5,
                help="plot width"),
    make_option(c("--height"), type="numeric", default=1.5,
                help="plot height")
)
usage <- "usage: %prog -i A.lst,B.lst -n TypeA,TypeB -o output.png\n"
parser <- OptionParser(option_list=option.list,
                       description="Draw Venn Diagram",
                       usage=usage)
opt <- parse_args(parser, args=args, positional_arguments=T)

inputFiles <- unlist(strsplit(opt$options$input, ",",
                       fixed = TRUE))
CatNames <- unlist(strsplit(opt$options$names, ",",
                     fixed = TRUE))
##print(inputFiles)
##str(inputFiles)
inputList <- list()
for(i in 1:length(inputFiles)){
    x <- read.table(inputFiles[i], header=F)
    inputList[[i]] <- x$V1
}
##str(input)
storedColor <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
if(length(inputFiles)==2){
    fillColor <- c(storedColor[1],storedColor[5])
}else if(length(inputFiles)==3){
    fillColor <- c(storedColor[1],storedColor[3],storedColor[5])
}else if(length(inputFiles)==4){
    fillColor <- c(storedColor[1],storedColor[2],storedColor[4],storedColor[5])
}else if(length(inputFiles)==5){
    fillColor <- storedColor
}else{
    write("Too many input lists", stderr())
}
p <- venn.diagram(inputList,
                  category.names = CatNames,
                  filename = NULL,
                  output = TRUE,
                  force.unique = TRUE,
                  ##imagetype = "png",
                  height = 720,
                  width = 720,
                  resolution = 300,
                  compression = "lzw",
                  lwd = 1,
                  main = opt$options$title,
                  main.cex = 0.6,
                  main.fontface = "bold", main.fontfamily = "sans",
                  fill = fillColor,
                  alpha = rep(0.65, length(inputList)),
                  cex = 0.4,
                  fontface = "bold",
                  fontfamily = "sans",
                  margin = 0.05,
                  cat.cex = 0.4, cat.fontface = "bold", cat.fontfamily = "sans")
ggsave(p, file=opt$options$output, width= opt$options$width, height = opt$options$height)