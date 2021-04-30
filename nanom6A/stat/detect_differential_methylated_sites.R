#! /usr/bin/env Rscript

library(tidyverse)
suppressPackageStartupMessages(library(DESeq2))

readData <- function(fileIn, sample){
    d <- read_tsv(fileIn,
                  col_names = c("gene", "chr", "pos", "motif", "ratio",
                                "methyl_count", "total_count")) %>%
        mutate(sample = sample)
    return(d)
}


sample1 <- readData("m6A_ratio.C08leaf_0h.tsv", "C08leaf_0h")
sample2 <- readData("m6A_ratio.C08leaf_1h.tsv", "C08leaf_1h")
sample3 <- readData("m6A_ratio.C08root_0h.tsv", "C08root_0h")
sample4 <- readData("m6A_ratio.C08root_1h.tsv", "C08root_1h")

ratio_tbl <- rbind(sample1, sample2,
                   sample3, sample4)
ratio_tbl <- ratio_tbl %>%
    mutate(tissue=str_sub(sample, 1, 7))


## Perform normalization
## Normalization will be performed within tissue
## `estimateSizeFactors` function of {DESeq2} was used
countNormalization <- function(ratio_tbl, selectedTissue){
    total_count_tbl <- ratio_tbl %>%
        mutate(id=paste(gene, chr, pos, sep=":")) %>%
        filter(tissue == selectedTissue) %>%
        select(id, sample, total_count) %>%
        pivot_wider(id_cols = id,
                    names_from = sample,
                    values_from = total_count) %>%
        column_to_rownames("id") %>%
        na.omit()

    colData <- data.frame(sample=c("sample1", "sample2"))
    rownames(colData) <- colnames(total_count_tbl)
    
    dds <- DESeqDataSetFromMatrix(countData = total_count_tbl,
                                  colData = colData,
                                  design = ~ sample)
    
    dds <- estimateSizeFactors(dds)

    total_count_tbl_normalized <- t(t(total_count_tbl)/sizeFactors(dds)) %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        pivot_longer(cols = !id,
                     names_to = "sample", values_to = "total_count")
    
    methyl_count_tbl <- ratio_tbl %>%
        mutate(id=paste(gene, chr, pos, sep=":")) %>%
        filter(tissue == selectedTissue) %>%
        select(id, sample, methyl_count) %>%
        pivot_wider(id_cols = id,
                    names_from = sample,
                    values_from = methyl_count) %>%
        column_to_rownames("id") %>%
        na.omit()
    
    methyl_count_tbl_normalized <- t(t(methyl_count_tbl)/sizeFactors(dds)) %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        pivot_longer(cols = !id,
                     names_to = "sample", values_to = "methyl_count")
    
    normalized_counts <- inner_join(methyl_count_tbl_normalized,
                                    total_count_tbl_normalized,
                                    by=c("id", "sample"))
    return(normalized_counts)
}

C08leaf_normalized <- countNormalization(ratio_tbl, "C08leaf")

C08root_normalized <- countNormalization(ratio_tbl, "C08root")


performFisherTest <- function(x){
    sample1_methylCount <- as.numeric(x[2])
    sample1_unmethylCount <- as.numeric(x[4])-as.numeric(x[2])
    sample2_methylCount <- as.numeric(x[3])
    sample2_unmethylCount <- as.numeric(x[5])-as.numeric(x[3])
    d <- data.frame(
        sample1=c(sample1_methylCount, sample1_unmethylCount),
        sample2=c(sample2_methylCount, sample2_unmethylCount)
    )
    pvalue <- fisher.test(d)$p.value
    return(pvalue)
}

C08leaf_ptable <- C08leaf_normalized %>%
    pivot_wider(id_cols = id,
                names_from = sample,
                values_from = c(methyl_count, total_count))

C08leaf_ptable$pvalue <- apply(C08leaf_ptable, 1, FUN=performFisherTest)

C08leaf_ptable  <- C08leaf_ptable %>%
    mutate(qvalue=p.adjust(pvalue, method = "BH")) %>%
    mutate(ratio_C08_leaf_0h=methyl_count_C08leaf_0h/total_count_C08leaf_0h,
           ratio_C08_leaf_1h=methyl_count_C08leaf_1h/total_count_C08leaf_1h)

C08root_ptable <- C08root_normalized %>%
    pivot_wider(id_cols = id,
                names_from = sample,
                values_from = c(methyl_count, total_count))

C08root_ptable$pvalue <- apply(C08root_ptable, 1, FUN=performFisherTest)

C08root_ptable <- C08root_ptable %>%
    mutate(qvalue=p.adjust(pvalue, method = "BH")) %>%
    mutate(ratio_C08_root_0h=methyl_count_C08root_0h/total_count_C08root_0h,
           ratio_C08_root_1h=methyl_count_C08root_1h/total_count_C08root_1h)


## Plot scatter plot for ratio
##
r1 <- ratio_tbl %>%
    select(gene, chr, pos, ratio, sample, tissue) %>%
    filter(tissue=="C08leaf") %>%
    pivot_wider(id_cols = c(gene, chr, pos, tissue),
                names_from = sample,
                values_from = ratio) %>%
    dplyr::rename(time_0h = C08leaf_0h,
                  time_1h = C08leaf_1h)

r2 <- ratio_tbl %>%
    select(gene, chr, pos, ratio, sample, tissue) %>%
    filter(tissue=="C08root") %>%
    pivot_wider(id_cols = c(gene, chr, pos, tissue),
                names_from = sample,
                values_from = ratio) %>%
    dplyr::rename(time_0h = C08root_0h,
                  time_1h = C08root_1h)

r_tbl <- rbind(r1, r2) %>% na.omit() %>%
    mutate(stat=case_when(
               time_1h - time_0h >= 0.1 ~ "up",
               time_1h - time_0h <=-0.1 ~ "down",
               TRUE ~ "unchanged"
           ),
           stat=factor(stat, levels = c("up", "unchanged", "down")))

count_tbl  <- r_tbl %>% group_by(tissue, stat) %>%
    summarize(count = n()) %>%
    ##mutate(y=rep(c(1,0.9,0.8), 2)) %>%
    mutate(x=0)
count_tbl$y <- rep(c(1,0.9,0.8), 2)

p <- ggplot(r_tbl,
            aes(x=time_0h,y=time_1h,
                color=stat))
p <- p + geom_point(alpha = 0.4)
p <- p + geom_text(data=count_tbl,
                   aes(x=x, y=y, label = paste0(stat, " : ", count)),
                   hjust = 0)
p <- p + xlab("Methylation level of m6A sites (0h)") + ylab("Methylation level of m6A sites (1h)")
p <- p + scale_color_manual(values = c("#ca0020", "darkgrey", "#0571b0"))
p <- p + facet_grid(~tissue)
p <- p + theme_bw(base_size=15)
p <- p + theme(legend.position = "none",
               axis.text = element_text(color = "black"),
               axis.ticks = element_line(color = "black"),
               panel.grid = element_blank(),
               strip.background = element_blank())

ggsave(p, filename = "m6A_ratio_scatter.pdf",
       width = 6, height = 3.5)

write_tsv(C08leaf_ptable, file = "C08leaf_dm6A.out.tsv")

write_tsv(C08root_ptable, file = "C08root_dm6A.out.tsv")

