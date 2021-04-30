#! /usr/bin/env Rscript

library(optparse)

## OptionParser
option_list <- list(
    make_option(c("-i", "--input"),
                help="Input GTF/GFF file"),
    make_option(c("-g", "--genome"),
                help="Path of genome fasta file"),
    make_option(c("-o", "--output"),
                help="Output File")
    )

opt <- parse_args(OptionParser(option_list=option_list,
                               description="Extract introns from GFF/GTF file, and detect canonical and non-canonical sites."))

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(tidyverse))

txdb <- makeTxDbFromGFF(opt$input)

## Generate Intron Granges object
introns <- intronicParts(txdb)

## Read reference sequences
refGenome <- readDNAStringSet(opt$genome)
intronSeq <- getSeq(refGenome, introns)

introns$donorSite = intronSeq %>% str_sub(1,2)
introns$accepterSite = intronSeq %>% str_sub(-2)

## Define canonical and non-canonical site
## Same criteria as STAR
## intron motif: 0: non-canonical;
## 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
introns$spliceSiteType = case_when(
    introns$donorSite == "GT" & introns$accepterSite == "AG" ~ "Canonical",
    introns$donorSite == "CT" & introns$accepterSite == "AC" ~ "Canonical",
    introns$donorSite == "GC" & introns$accepterSite == "AG" ~ "Canonical",
    introns$donorSite == "CT" & introns$accepterSite == "GC" ~ "Canonical",
    introns$donorSite == "AT" & introns$accepterSite == "AC" ~ "Canonical",
    introns$donorSite == "GT" & introns$accepterSite == "AT" ~ "Canonical",
    TRUE ~ "Non-Canonical"
)

introns$tx_name <- paste(introns$tx_name, collapse = ",")

introns <- as.data.frame(introns)
##head(introns)
write.table(introns,
            file = opt$output,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)
