#! /usr/bin/env Rscript

library(tidyverse)
library(QuantifyPolyA)

bedFiles <- dir(pattern = "^(Leaf|Root)_[01]h_R[123].PAS.bed")

QpolyA <- Load.PolyA(files= bedFiles)

QpolyA <- Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 12)

QpolyA <- Annotate.PolyA(QpolyA,
                         gff='/home/ubuntu/salinity_suppl_analysis/analysis/transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf',
                         seq.levels = NA)


QpolyA <- Filter.PolyA(QpolyA, min_count = 5, min_sample = 1)

## write filtered result to bed file for JBrowse
full_polyA_sites <- QpolyA@polyA %>%
    rownames_to_column("polyA_id") %>%
    mutate(siteID = paste(polyA_id, gene_id, type, score, sep="::")) %>%
    dplyr::select(seqnames, start, end, siteID, strand) %>%
    mutate(bedScore=0) %>%
    dplyr::relocate(bedScore, .before = strand)

full_polyA_sites %>% write_tsv("QuantifyPolyA_polyA_sites.bed",
                               col_names = FALSE)

sampleInfo <- data.frame(
    condition = str_replace(QpolyA@sample_names,
                            "_R[123].PAS",""),
    replicate = str_extract(QpolyA@sample_names,
                            "R[123]")
)
rownames(sampleInfo) <- QpolyA@sample_names

## Generate PCA
res <- DESeq2.PolyA(QpolyA, sampleInfo)

p <- res$PCA.Plot + scale_color_brewer(palette = "Dark2")+
    theme(panel.grid = element_blank())
ggplot2::ggsave(plot = p, filename = "polyA_PCA.pdf",
                width = 5,
                height = 4)


## Extract upstream 50 nt sequences for motif discovery
QpolyA@polyA %>%
    dplyr::select(seqnames, strand, center) %>%
    as_tibble() %>%
    mutate(
        start = case_when(
            `strand` == "+" ~ center-49,
            `strand` == "-" ~ center+1
        ),
        end = case_when(
            `strand` == "+" ~ center+0,
            `strand` == "-" ~ center+50
        )
    ) %>%
    mutate(id=rownames(QpolyA@polyA), bedScore = 0) %>%
    dplyr::select(seqnames, start, end, id, bedScore, strand) %>%
    write_tsv(file = "polyA_upstream50nt.bed", col_names = FALSE)

## Extract upstream 200 nt for motif occurrence illustration
QpolyA@polyA %>%
    dplyr::select(seqnames, strand, center) %>%
    as_tibble() %>%
    mutate(
        start = case_when(
            `strand` == "+" ~ center-199,
            `strand` == "-" ~ center+1
        ),
        end = case_when(
            `strand` == "+" ~ center+0,
            `strand` == "-" ~ center+200
        )
    ) %>%
    mutate(id=rownames(QpolyA@polyA), bedScore = 0) %>%
    dplyr::select(seqnames, start, end, id, bedScore, strand) %>%
    write_tsv(file = "polyA_upstream200nt.bed", col_names = FALSE)


## Conduct dynamic test
gene.APA.salt_root <- Quantify.GeneAPA(
    QpolyA,
    sampleInfo,
    contrast=c("condition","Root_0h","Root_1h")
)


gene.APA.salt_leaf <- Quantify.GeneAPA(
    QpolyA,
    sampleInfo,
    contrast=c("condition", "Leaf_0h", "Leaf_1h")
)

gene.APA.root_vs_leaf_0h <- Quantify.GeneAPA(
    QpolyA,
    sampleInfo,
    contrast=c("condition","Leaf_0h","Root_0h")
)

gene.APA.root_vs_leaf_1h <- Quantify.GeneAPA(
    QpolyA,
    sampleInfo,
    contrast=c("condition","Leaf_1h","Root_1h")
)

no.na <- function(x) ifelse(is.na(x), 1, x)

gene.APA.salt_root$p.value <- no.na(gene.APA.salt_root$p.value)

gene.APA.salt_root$padj <- p.adjust(gene.APA.salt_root$p.value, method = "BH")

gene.APA.salt_leaf$p.value <- no.na(gene.APA.salt_leaf$p.value)

gene.APA.salt_leaf$padj <- p.adjust(gene.APA.salt_leaf$p.value, method = "BH")

gene.APA.root_vs_leaf_0h$p.value <- no.na(gene.APA.root_vs_leaf_0h$p.value)

gene.APA.root_vs_leaf_0h$padj <- p.adjust(gene.APA.root_vs_leaf_0h$p.value, method = "BH")

gene.APA.root_vs_leaf_1h$p.value <- no.na(gene.APA.root_vs_leaf_1h$p.value)

gene.APA.root_vs_leaf_1h$padj <- p.adjust(gene.APA.root_vs_leaf_1h$p.value, method = "BH")

write_tsv(gene.APA.salt_leaf, file = "salt_leaf_1hvs0h.gene_APA.tsv")
write_tsv(gene.APA.salt_root, file = "salt_root_1hvs0h.gene_APA.tsv")
write_tsv(gene.APA.root_vs_leaf_0h,
          file = "salt_root_vs_leaf_0h.gene_APA.tsv")
write_tsv(gene.APA.root_vs_leaf_1h,
          file = "salt_root_vs_leaf_1h.gene_APA.tsv")


## Count lengthening and shortening events
gene.APA.salt_root  <- gene.APA.salt_root %>%
    mutate(stat=case_when(
               r>0 ~ "lengthening",
               r<0 ~ "shortening",
               TRUE ~ "nochange")
           )

salt_root_stats <- gene.APA.salt_root %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    group_by(stat) %>% summarise(count=n())

salt_root_lengthening <- gene.APA.salt_root %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    dplyr::filter(stat == "lengthening") %>%
    dplyr::select(gene_id)

salt_root_shortening <- gene.APA.salt_root %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    dplyr::filter(stat == "shortening") %>%
    dplyr::select(gene_id)

gene.APA.salt_leaf <- gene.APA.salt_leaf %>%
    mutate(stat=case_when(
               r>0 ~ "lengthening",
               r<0 ~ "shortening",
               TRUE ~ "nochange")
           )

salt_leaf_stats <- gene.APA.salt_leaf %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    group_by(stat) %>% summarise(count=n())

salt_leaf_lengthening <- gene.APA.salt_leaf %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    dplyr::filter(stat == "lengthening") %>%
    dplyr::select(gene_id)

salt_leaf_shortening <- gene.APA.salt_leaf %>%
    dplyr::filter(!is.na(r), !is.nan(r)) %>%
    dplyr::filter(stat == "shortening") %>%
    dplyr::select(gene_id)

write_tsv(salt_leaf_lengthening,
          "salt_leaf_1hvs0h.lengthening.lst",
          col_names = FALSE)
write_tsv(salt_leaf_shortening,
          "salt_leaf_1hvs0h.shortening.lst",
          col_names = FALSE)
write_tsv(salt_root_lengthening,
          "salt_root_1hvs0h.lengthening.lst",
          col_names = FALSE)
write_tsv(salt_root_shortening,
          "salt_root_1hvs0h.shortening.lst",
          col_names = FALSE)



anno = import("/home/ubuntu/salinity_suppl_analysis/analysis/transdecoder/stringtie.count_5.withCDS.strandCorrected.gtf")

gene_id = "XLOC_011380"

dp = Visualize.PolyA(QpolyA,
                     anno, gene_id, coord.lim=c(35453250,35460010))

dp + scale_color_aaas() + scale_fill_aaas()

ggsave(plot = dp, filename = "example.pdf")

save.image("QuantifyPolyA.RData")
