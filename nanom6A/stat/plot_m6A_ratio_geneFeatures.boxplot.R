#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))

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

txdb <- makeTxDbFromGFF("../../transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf")

##selectedTrans <- selectLongestTranscript(txdb, TRUE)
selectedTrans <- selectLongestTranscript(txdb, FALSE)

## Generate reference GRs from txdb
## fiveUTR
fiveUTR <- fiveUTRsByTranscript(txdb, use.names = T)
fiveUTR_txList <- intersect(names(fiveUTR), selectedTrans)
fiveUTR <- fiveUTR[fiveUTR_txList]
## threeUTR
threeUTR <- threeUTRsByTranscript(txdb, use.names = T)
threeUTR_txList <- intersect(names(threeUTR), selectedTrans)
threeUTR <- threeUTR[threeUTR_txList]
## CDS
cds <- cdsBy(txdb, by="tx", use.name = T)
cds_txList <- intersect(names(cds), selectedTrans)
cds <- cds[cds_txList]

findSite_in_refGR <- function(query_gr, subject_gr){
    query_gr[queryHits(GenomicRanges::findOverlaps(m6A_sites_gr,
                                                   subject_gr))]
}

## Generate data for C08leaf
d <- read_tsv("C08leaf_dm6A.out.tsv") %>%
    dplyr::select(id, ratio_C08_leaf_0h, ratio_C08_leaf_1h) %>%
    separate(id, into = c("geneID", "chr", "pos"),
             remove=FALSE, sep = ":")

m6A_sites_gr <- GRanges(
    seqnames = d$chr,
    ranges = IRanges(d$pos),
    id = d$id,
    ratio_C08_leaf_0h = d$ratio_C08_leaf_0h,
    ratio_C08_leaf_1h = d$ratio_C08_leaf_1h
)

fiveUTR_m6A <- findSite_in_refGR(m6A_sites_gr, fiveUTR)
threeUTR_m6A <- findSite_in_refGR(m6A_sites_gr, threeUTR)
cds_m6A <- findSite_in_refGR(m6A_sites_gr, cds)

## Init tibble
m6A_ratio_Long <- tibble(
    id = character(0),
    chr = character(0),
    pos = numeric(0),
    ratio = numeric(0),
    sample= character(0),
    feature= character(0)
)

## Populate data
m6A_ratio_Long  <- m6A_ratio_Long %>%
    bind_rows(
        tibble(
            id = elementMetadata(fiveUTR_m6A)[,1],
            chr = as.vector(seqnames(fiveUTR_m6A)),
            pos = start(fiveUTR_m6A),
            ratio = elementMetadata(fiveUTR_m6A)[,2],
            sample= "C08leaf_0h",
            feature= "5' UTR"
        ),
        tibble(
            id = elementMetadata(fiveUTR_m6A)[,1],
            chr = as.vector(seqnames(fiveUTR_m6A)),
            pos = start(fiveUTR_m6A),
            ratio = elementMetadata(fiveUTR_m6A)[,3],
            sample= "C08leaf_1h",
            feature= "5' UTR"
        ),
        tibble(
            id = elementMetadata(threeUTR_m6A)[,1],
            chr = as.vector(seqnames(threeUTR_m6A)),
            pos = start(threeUTR_m6A),
            ratio = elementMetadata(threeUTR_m6A)[,2],
            sample= "C08leaf_0h",
            feature= "3' UTR"
        ),
        tibble(
            id = elementMetadata(threeUTR_m6A)[,1],
            chr = as.vector(seqnames(threeUTR_m6A)),
            pos = start(threeUTR_m6A),
            ratio = elementMetadata(threeUTR_m6A)[,3],
            sample= "C08leaf_1h",
            feature= "3' UTR"
        ),
        tibble(
            id = elementMetadata(cds_m6A)[,1],
            chr = as.vector(seqnames(cds_m6A)),
            pos = start(cds_m6A),
            ratio = elementMetadata(cds_m6A)[,2],
            sample= "C08leaf_0h",
            feature= "CDS"
        ),
        tibble(
            id = elementMetadata(cds_m6A)[,1],
            chr = as.vector(seqnames(cds_m6A)),
            pos = start(cds_m6A),
            ratio = elementMetadata(cds_m6A)[,3],
            sample= "C08leaf_1h",
            feature= "CDS"
        )
    )

## Generate data for C08root
d <- read_tsv("C08root_dm6A.out.tsv") %>%
    dplyr::select(id, ratio_C08_root_0h, ratio_C08_root_1h) %>%
    separate(id, into = c("geneID", "chr", "pos"),
             remove=FALSE, sep = ":")

m6A_sites_gr <- GRanges(
    seqnames = d$chr,
    ranges = IRanges(d$pos),
    id = d$id,
    ratio_C08_root_0h = d$ratio_C08_root_0h,
    ratio_C08_root_1h = d$ratio_C08_root_1h
)

fiveUTR_m6A <- findSite_in_refGR(m6A_sites_gr, fiveUTR)
threeUTR_m6A <- findSite_in_refGR(m6A_sites_gr, threeUTR)
cds_m6A <- findSite_in_refGR(m6A_sites_gr, cds)

## Populate data
m6A_ratio_Long  <- m6A_ratio_Long %>%
    bind_rows(
        tibble(
            id = elementMetadata(fiveUTR_m6A)[,1],
            chr = as.vector(seqnames(fiveUTR_m6A)),
            pos = start(fiveUTR_m6A),
            ratio = elementMetadata(fiveUTR_m6A)[,2],
            sample= "C08root_0h",
            feature= "5' UTR"
        ),
        tibble(
            id = elementMetadata(fiveUTR_m6A)[,1],
            chr = as.vector(seqnames(fiveUTR_m6A)),
            pos = start(fiveUTR_m6A),
            ratio = elementMetadata(fiveUTR_m6A)[,3],
            sample= "C08root_1h",
            feature= "5' UTR"
        ),
        tibble(
            id = elementMetadata(threeUTR_m6A)[,1],
            chr = as.vector(seqnames(threeUTR_m6A)),
            pos = start(threeUTR_m6A),
            ratio = elementMetadata(threeUTR_m6A)[,2],
            sample= "C08root_0h",
            feature= "3' UTR"
        ),
        tibble(
            id = elementMetadata(threeUTR_m6A)[,1],
            chr = as.vector(seqnames(threeUTR_m6A)),
            pos = start(threeUTR_m6A),
            ratio = elementMetadata(threeUTR_m6A)[,3],
            sample= "C08root_1h",
            feature= "3' UTR"
        ),
        tibble(
            id = elementMetadata(cds_m6A)[,1],
            chr = as.vector(seqnames(cds_m6A)),
            pos = start(cds_m6A),
            ratio = elementMetadata(cds_m6A)[,2],
            sample= "C08root_0h",
            feature= "CDS"
        ),
        tibble(
            id = elementMetadata(cds_m6A)[,1],
            chr = as.vector(seqnames(cds_m6A)),
            pos = start(cds_m6A),
            ratio = elementMetadata(cds_m6A)[,3],
            sample= "C08root_1h",
            feature= "CDS"
        )
    )



m6A_ratio_Wide <- m6A_ratio_Long %>% unique() %>%
    pivot_wider(id_cols = c(id, chr, pos, feature),
                names_from = sample,
                values_from = ratio)

m6A_deltaRatio_Long <- m6A_ratio_Wide %>%
    mutate(delta_C08leaf=C08leaf_1h-C08leaf_0h,
           delta_C08root=C08root_1h-C08root_0h) %>%
    dplyr::select(chr,pos,feature,delta_C08leaf,delta_C08root) %>%
    pivot_longer(4:5, names_to = "sample",
                 values_to = "deltaRatio",
                 names_prefix="delta_") %>%
    na.omit()

m6A_deltaRatio_Long <- m6A_deltaRatio_Long %>%
    mutate(
        stat=case_when(
            deltaRatio>=0.1 ~ "up",
            deltaRatio<=-0.1 ~ "down",
            TRUE ~ "unchanged"
        )
    )

p <- ggplot(m6A_deltaRatio_Long,
            aes(x=feature, y=deltaRatio, fill=feature)) +
    geom_violin(alpha=0.8) +
    geom_boxplot(width=0.2, outlier.shape = NA, alpha = 0.8) +
    facet_grid(stat~sample, scales = "free")

p <- p + theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(color="black"),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust =1))

ggsave(filename = "m6A_ratio_GenomicFeatures.pdf",
       width = 4, height = 3.5)

