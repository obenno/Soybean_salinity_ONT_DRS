#! /usr/bin/env Rscript

library(optparse)

option_list <- list(
    make_option(c("-q", "--query"), type = "character",
                help = "Query GTF file"),
    make_option(c("-s", "--subject"), type = "character",
                help = "Subject GTF file"),
    make_option(c("-t", "--threads"), type="integer", default = 8,
                help = "Number of threads to use"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output filename")
)

parser <- OptionParser(option_list = option_list,
                       description = "Detect fusion transcripts by comparing new gtf to reference gtf.")

opt <- parse_args(parser)

if(length(commandArgs(trailingOnly = T))==0){
    print_help(parser)
    quit("no")
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(doParallel))

##numCores <- detectCores()
##registerDoParallel(numCores)
registerDoParallel(opt$threads)

txdb1 <- makeTxDbFromGFF(opt$query)

query_exonsList <- exonsBy(txdb1, by="tx", use.names = T)

txdb2 <- makeTxDbFromGFF(opt$subject)

subject_exons <- exons(txdb2, columns = c("exon_id", "gene_id"))

##fusion_tx <- tibble(id=character(0), geneID=character(0))

fusion_tx <- foreach(i=1:length(query_exonsList), .combine=rbind) %dopar% {
    overlap_result <- findOverlaps(query_exonsList[[i]], subject_exons,
                                   type = "equal")
    gene_hits <- mcols(subject_exons[subjectHits(overlap_result)])$gene_id %>%
                                                     unlist() %>%
                                                     unique()
    if(length(gene_hits)>1){
        fusionID <- names(query_exonsList)[i]
        geneIDs <- paste0(gene_hits, collapse = ",")
        out <- tibble_row(id=fusionID, geneID=geneIDs)
        ## fusion_tx <- fusion_tx %>%
        ##     bind_rows(out)
    }
}

write_tsv(fusion_tx, file = opt$output)

