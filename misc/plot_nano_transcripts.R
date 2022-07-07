#! /usr/bin/env Rscript

library(optparse)

## OptionParser
option_list <- list(
    make_option(c("-i", "--input"),
                help="Input GTF/GFF file as reference"),
    make_option(c("-c", "--coordinate"), type = "character",
                help="Genomic coordinate to plot (Chr01:1000-2000)"),
    make_option(c("--highlight"), type = "character", default=NULL,
                help="Highlight genomic region (Chr01:1000-2000)"),
    make_option(c("--anno"), default=NULL,
                help="Annotation data (GTF/BED)"),
    make_option(c("--align"), default=NULL,
                help="Alignment data (bam file)"),
    make_option(c("-o", "--output"),
                help="Output File"),
    make_option(c("--width"), type="double",
                help="Output file width"),
    make_option(c("--height"), type="double",
                help="Output file height"),
    make_option(c("--sizes"), default=NULL,
                help="track heights, comma separated")
)

parser <- OptionParser(option_list=option_list,
                       description="Plot gene regions with Gviz")

opt <- parse_args(parser)

if(length(commandArgs(trailing=TRUE))==0){
    print_help(parser)
    quit("no")
}

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(tidyverse))

options(ucscChromosomeNames=FALSE)

txdb <- makeTxDbFromGFF(opt$input)

coordinate <- str_split(opt$coordinate, ":",
                        simplify = TRUE)
chr <- coordinate[1]
region <- str_split(coordinate[2], "-",
                    simplify = TRUE)
start <- as.numeric(region[1])
end <- as.numeric(region[2])

## Init plotTrack
allTracks <- list()
## Axis track
axisTrack <- GenomeAxisTrack(col="black",
                             fontsize=6,
                             fontcolor="black",
                             lwd=1,
                             distFromAxis=0.5)
allTracks <- c(allTracks, axisTrack)
## Reference gene track
refTrack <- GeneRegionTrack(txdb,
                            chromosome = chr,
                            start=start, end=end,
                            shape=c("smallArrow", "box", "arrow", "ellipse", "fixedArrow"),
                            col = NULL,
                            col.line=NULL,
                            fill = "#1b9e77",
                            arrowHeadMaxWidth=40,
                            arrowHeadWidth=20)

allTracks <- c(allTracks, refTrack)
message("Ref track created.")
## Annotation track(s)
if(!is.null(opt$anno)){
    annoFiles <- str_split(opt$anno, ",",
                           simplify = TRUE)
    annoTracks <- list()
    for(i in 1:length(annoFiles)){
        if(str_detect(annoFiles[i], "bed$")){
            annoTracks[i] <- AnnotationTrack(annoFiles[i],
                                             chromosome = chr,
                                             start=start, end=end)
        }else if(str_detect(annoFiles[i], "gtf$")){
            d <- makeTxDbFromGFF(annoFiles[i])
            annoTracks[i] <- GeneRegionTrack(d,
                                             chromosome = chr,
                                             start=start, end=end,
                                             col=NULL,
                                             col.line=NULL,
                                             fill = "#7570b3")
        }else{
            stop("Format not supported.")
        }
    }
}else{
    annoTracks <- NULL
}
allTracks <- c(allTracks, annoTracks)
message("Anno tracks created.")

## Mapping alignment track
if(!is.null(opt$align)){
    alignFiles <- str_split(opt$align, ",", simplify = TRUE)
    alignTracks <- list()
    for(i in 1:length(alignFiles)){
        alignTracks[i] <- AlignmentsTrack(alignFiles[i],
                                         chromosome = chr,
                                         start=start,
                                         end=end,
                                         isPaired = FALSE,
                                         alpha.reads = 0.8,
                                         lwd.reads = 0.5,
                                         col.reads="darkgrey",
                                         type = "pileup")
    }
}else{
    alignTracks <- NULL
}

allTracks <- c(allTracks, alignTracks)
message("Bam track created.")

message("In total ", length(allTracks), " tracks.")
sizes <- str_split(opt$sizes, ",", simplify = TRUE) %>% as.numeric()

if(!is.null(opt$highlight)){
    message("Hightlight region: ", opt$highlight)
    ht_coordinate <- str_split(opt$highlight, ":",
                               simplify =  TRUE)
    ht_chr <- ht_coordinate[1]
    ht_region <- str_split(ht_coordinate[2], "-",
                           simplify =  TRUE)
    ht_start <- as.numeric(ht_region[1])
    ht_end <- as.numeric(ht_region[2])
    ht <- HighlightTrack(allTracks[2:length(allTracks)],
                         start = ht_start, width = ht_end - ht_start +1,
                         chromosome = ht_chr,
                         col = "#d95f02",
                         fill = "#fed4b4")
    allTracks <- c(allTracks[1], ht)
}

## ## Set plot scheme
## scheme <- getScheme()
## scheme$GeneRegionTrack$fill <- "salmon"
## scheme$GeneRegionTrack$col <- NULL
## addScheme(scheme, "myScheme")
## options(Gviz.scheme = "myScheme")

pdf(opt$output, width = opt$width,
    height = opt$height)
plotTracks(allTracks,
           chromosome = chr,
           from = start,
           to=end,
           sizes=sizes,
           showTitle = FALSE)

dev.off()

