#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-i", "--inputGFF"), default=NULL,
                help="Print extra output [default]"),
    make_option(c("-s", "--signal"), default=NULL,
                help="Features to be summrized, bed format, comma separated list"),
    make_option(c("-n", "--names"), default=NULL, dest="sample",
                help="Sample names, comma separated list"),
    make_option(c("-l", "--level"), default="gene",
                help="Statistic on gene or transcript level [%default]"),
    make_option(c("-o", "--output"), default=NULL,
                help="Output plot filename"),
    make_option(c("-t", "--threads"), default=NULL, type="integer",
                help="Threads to be used [%default]")
)

parser <- OptionParser(option_list=option_list)

opt <- parse_args(parser)

if(length(commandArgs(trailing=TRUE))==0){
    print_help(parser)
    quit("no")
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(doParallel))

if(is.null(opt$threads)){
    numCores <- detectCores()
}else{
    numCores <- opt$threads
}

message("Used parallel with ", numCores, " thread(s).")
registerDoParallel(numCores)


##signal <- import(opt$signal)

## fiveUTR <- fiveUTRsByTranscript(txdb, use.names = T)
## threeUTR <- threeUTRsByTranscript(txdb, use.names = T)
## cds <- cdsBy(txdb, by="tx", use.name = T)
##
## fiveUTR_coverage <- coverageByTranscript(signal, fiveUTR)

selectLongestTranscript <- function(txdb, switch){
    d <- transcriptLengths(txdb) %>% arrange(gene_id, desc(tx_len))
    if(switch){
        gene <- c()
        selectedTrans <- c()
        for(i in 1:nrow(d)){
            if(!(d[i,"gene_id"] %in% gene)){
                selectedTrans <- c(selectedTrans, d[i,"tx_name"])
                gene <- c(gene, d[i,"gene_id"])
            }
        }
    }else{
        selectedTrans <- d$tx_name
    }
    return(selectedTrans)
}


if(opt$level=="gene"){
    ## Select longest transcript for each gene
    ## as represented ones
    switch <- TRUE
}else{
    switch <- FALSE
}



bin_and_count <- function(Rle, binNum){
    binned_Rle = split(as.numeric(Rle), cut(seq_along(as.numeric(Rle)), binNum))
    featureCounts <- map_df(binned_Rle, function(x) sum(x)) ##%>% as.numeric()
    names(featureCounts)  <- paste0("bin_", c(1:binNum))
    return(featureCounts)
}



GenerateSampleCount <- function(txdb, selectedTrans, signal,
                                fiveUTR_binNum, threeUTR_binNum,
                                cds_binNum){
    ## Start fiveUTR

    fiveUTR <- fiveUTRsByTranscript(txdb, use.names = T)
    fiveUTR_txList <- intersect(names(fiveUTR), selectedTrans)
    fiveUTR <- fiveUTR[fiveUTR_txList]

    fiveUTR_coverage <- coverageByTranscript(signal, fiveUTR)

    k=foreach (i=1:length(fiveUTR_coverage)) %dopar% {
        bin_and_count(fiveUTR_coverage[[i]], fiveUTR_binNum)
    }

    fiveUTR_count_tbl <- bind_rows(k) %>%
        mutate(id=names(fiveUTR_coverage)) %>%
        pivot_longer(cols = -id,
                     names_to = "binName",
                     values_to = "count") %>%
        mutate(type="fiveUTR") %>%
        group_by(type, binName) %>%
        summarize(count=sum(count))

    ## Start threeUTR

    threeUTR <- threeUTRsByTranscript(txdb, use.names = T)
    threeUTR_txList <- intersect(names(threeUTR), selectedTrans)
    threeUTR <- threeUTR[threeUTR_txList]

    threeUTR_coverage <- coverageByTranscript(signal, threeUTR)

    k=foreach (i=1:length(threeUTR_coverage)) %dopar% {
        bin_and_count(threeUTR_coverage[[i]], threeUTR_binNum)
    }

    threeUTR_count_tbl <- bind_rows(k) %>%
        mutate(id=names(threeUTR_coverage)) %>%
        pivot_longer(cols = -id,
                     names_to = "binName",
                     values_to = "count") %>%
        mutate(type="threeUTR") %>%
        group_by(type, binName) %>%
        summarize(count=sum(count))

    ## Start CDS

    cds <- cdsBy(txdb, by="tx", use.name = T)
    cds_txList <- intersect(names(cds), selectedTrans)
    cds <- cds[cds_txList]

    cds_coverage <- coverageByTranscript(signal, cds)

    k=foreach (i=1:length(cds_coverage)) %dopar% {
        bin_and_count(cds_coverage[[i]], cds_binNum)
    }

    cds_count_tbl <- bind_rows(k) %>%
        mutate(id=names(cds_coverage)) %>%
        pivot_longer(cols = -id,
                     names_to = "binName",
                     values_to = "count") %>%
        mutate(type="CDS") %>%
        group_by(type, binName) %>%
        summarize(count=sum(count))


    count_tbl <- rbind(fiveUTR_count_tbl,
                       threeUTR_count_tbl,
                       cds_count_tbl) %>%
        mutate(featureBin=paste(type, binName, sep=":"))

    bin_order <- paste0("bin_", c(1:fiveUTR_binNum, 1:cds_binNum, 1:threeUTR_binNum))

    x_order <- paste(c(rep("fiveUTR", fiveUTR_binNum), rep("CDS", cds_binNum), rep("threeUTR",threeUTR_binNum)), bin_order, sep=":")

    count_tbl  <- count_tbl %>%
        mutate(featureBin=factor(featureBin,
                                 levels = x_order))
    return(count_tbl)
}

txdb <- makeTxDbFromGFF(opt$input)
signalFiles <- unlist(str_split(opt$signal, ","))
samples <- unlist(str_split(opt$sample, ","))

selectedTrans <- selectLongestTranscript(txdb, switch)

fiveUTR_binNum <- 20
threeUTR_binNum <- 20
cds_binNum <- 60

count_tbl <- list()
for(i in 1:length(signalFiles)){
    signal <- import(signalFiles[i])
    count_tbl[[i]] <- GenerateSampleCount(txdb, selectedTrans,
                                          signal,
                                          fiveUTR_binNum, threeUTR_binNum,
                                          cds_binNum)
    count_tbl[[i]] <- count_tbl[[i]] %>%
        mutate(sample=samples[i])
}

count_tbl <- bind_rows(count_tbl)

count_tbl <- count_tbl %>%
    dplyr::select(featureBin, sample, count) %>%
    tidyr::uncount(count)

if(length(samples)>2){
    legend_y <- 0.8
}else{
    legend_y <- 0.85
}

p <- ggplot(count_tbl,
            aes(x=featureBin,## y=count,
                group=sample,
                fill=sample))
p <- p + geom_density(alpha = 0.6)
p <- p + xlab("") + ylab("m6A sites density")
##p <- p + geom_bar(stat = "identity", alpha=0.6,
##                  position = "identity")
p <- p + scale_fill_brewer(palette = "Dark2")
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank(),
               legend.position = c(0.1, legend_y),
               legend.title = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank())

p <- p + geom_vline(xintercept= "CDS:bin_1", linetype = "dashed")

p <- p + geom_vline(xintercept= paste0("CDS:bin_", cds_binNum),
                    linetype = "dashed")

p <- p + annotate("rect", fill="darkgrey", size=0.5,
                  color="black",
                  xmin="CDS:bin_1",
                  xmax = paste0("CDS:bin_", cds_binNum),
                  ymin = -0.001, ymax = -0.0001)

p <- p + annotate("segment", x="fiveUTR:bin_1",
                  xend = "CDS:bin_1",
                  y=-0.0005, yend=-0.0005, size=1.5)

p <- p + annotate("segment", x=paste0("CDS:bin_", cds_binNum),
                  xend = paste0("threeUTR:bin_", threeUTR_binNum),
                  y=-0.0005, yend=-0.0005, size=1.5)

p <- p + annotate("text", x=paste0("fiveUTR:bin_", round(fiveUTR_binNum/2)),
                  y= 0.001, label="5' UTR",
                  fontface = "bold")

p <- p + annotate("text", x=paste0("threeUTR:bin_", round(threeUTR_binNum/2)),
                  y= 0.001, label="3' UTR",
                  fontface = "bold")

p <- p + annotate("text", x=paste0("CDS:bin_", round(cds_binNum/2)),
                  y= 0.001, label="CDS",
                  fontface = "bold")

ggsave(filename = opt$output,
       height = 3.5,
       width = 7)

