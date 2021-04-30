#! /usr/bin/env Rscript

library(tidyverse)

C08leaf_0h_polyA <- read_tsv("transcript_polyA_length.C08leaf_0h.tsv",
                          col_names = c('transID', 'readID',
                                        'polyA_length')) %>%
    mutate(sample = 'C08leaf_0h')


C08leaf_1h_polyA <- read_tsv("transcript_polyA_length.C08leaf_1h.tsv",
                             col_names = c('transID', 'readID',
                                           'polyA_length')) %>%
    mutate(sample = 'C08leaf_1h')

C08root_0h_polyA <- read_tsv("transcript_polyA_length.C08root_0h.tsv",
                             col_names = c('transID', 'readID',
                                           'polyA_length')) %>%
    mutate(sample = 'C08root_0h')

C08root_1h_polyA <- read_tsv("transcript_polyA_length.C08root_1h.tsv",
                             col_names = c('transID', 'readID',
                                           'polyA_length')) %>%
    mutate(sample = 'C08root_1h')

polyA_profile <- rbind(C08leaf_0h_polyA, C08leaf_1h_polyA,
                       C08root_0h_polyA, C08root_1h_polyA)

calc_polyA_median <- function(polyA_string){
    reads_polyA_len <- unlist(str_split(polyA_string, ",")) %>%
        map_dbl(.f=as.numeric)
    return(median(reads_polyA_len))
}

polyA_profile <- polyA_profile %>%
    mutate(polyA_median = map_dbl(.x=polyA_length,
                                  .f=calc_polyA_median)) %>%
    mutate(tissue = str_sub(sample, 1,7))

## Disabled, not necessary any more

## Stats for transcripts reads count
## trans_reads <- polyA_profile %>%
##     mutate(readsList = map(readID, .f=function(x) unlist(str_split(x, ","))),
##            readsCount = map(readsList, length)) %>%
##     select(transID, sample, readsCount) %>% unnest(readsCount)
## 
## p1 <- ggplot(trans_reads, aes(x=readsCount))
## p1 <- p1 + geom_histogram(binwidth = 1)
## p1 <- p1 + facet_grid(~sample)
## p1 <- p1 + xlim(4,50)
## ggsave(p1, filename = "trans_polyAreads_histogram.pdf",
##        width = 6,
##        height = 3)

## polyA_profile %>%
##     filter(tissue == "C08root") %>%
##     kruskal.test(polyA_median ~ sample, data = .data)


## Mann-Whitney test between treatment for both tissue
x <- polyA_profile %>%
    filter(sample == "C08root_1h") %>% pull(polyA_median)
y <- polyA_profile %>%
    filter(sample == "C08root_0h") %>% pull(polyA_median)
C08root_pvalue <- wilcox.test(x, y)$p.value %>%
                                  signif(digits=3)

## polyA_profile %>%
##     filter(tissue == "C08leaf") %>%
##     as.data.frame() %>%
##     kruskal.test(polyA_median ~ sample)

x <- polyA_profile %>%
    filter(sample == "C08leaf_1h") %>% pull(polyA_median)
y <- polyA_profile %>%
    filter(sample == "C08leaf_0h") %>% pull(polyA_median)
C08leaf_pvalue <- wilcox.test(x, y)$p.value %>%
                                  signif(digits=3)


p_tbl <- tibble(tissue = c("C08leaf", "C08root"),
                pvalue = c(C08leaf_pvalue, C08root_pvalue))

count_tbl <- polyA_profile %>%
    group_by(sample, tissue) %>%
    summarise(count=n())

p <- ggplot(polyA_profile,
            aes(x=sample, y=polyA_median,
                fill=sample))

## p <- p + stat_bin(position = "identity",
##                   alpha = 0.6,
##                   bins = 100)
p <- p + geom_violin(alpha = 0.8,
                     trim = TRUE)
p <- p + geom_boxplot(width = 0.2,
                      alpha = 0.8,
                      outlier.shape = NA)

p <- p + geom_text(data = count_tbl,
                   aes(x=sample, y=380,
                       label = paste0("n = ", count)),
                   size = 3)

p <- p + geom_text(data = p_tbl,
                   aes(x=1.5, y=300, fill = NULL,
                       label = paste0("P < ", pvalue)),
                   size = 3)

p <- p + geom_segment(aes(x = 1.1, y = 280,
                          xend = 1.9, yend = 280))

p <- p + facet_grid(~tissue, scales = "free")
##p <- p + ylim(0,300)
p <- p + xlab("") + ylab("Transcripts polyA length (nt)")
p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               strip.background = element_blank(),
               axis.text = element_text(color = "black"))

ggsave(filename = "polyA_histogram.pdf",
       width = 5,
       height = 3)

## Format profile and combine data from
## two condition (0h, 1h)
## keep the transcripts detected in both samples

compare_polyA_between_samples <- function(data_longFormat,
                                          sample1, sample2){
    out <- data_longFormat %>%
        filter(sample == sample1 | sample == sample2) %>%
        select(transID, polyA_length, sample) %>%
        pivot_wider(names_from = sample, values_from = polyA_length) %>%
        na.omit() %>%
        mutate(sample1_data = str_split(.data[[sample1]], ","),
               sample2_data = str_split(.data[[sample2]], ",")) %>%
        mutate(sample1_data = map(sample1_data, as.numeric),
               sample2_data = map(sample2_data, as.numeric)) %>%
        mutate(sample1_median = map_dbl(sample1_data, median)) %>%
        mutate(sample2_median = map_dbl(sample2_data, median)) %>%
        mutate(pvalue = map2_dbl(.x = sample1_data,
                                 .y = sample2_data,
                                 .f = function(x1, x2) wilcox.test(x1, x2, exact = FALSE)$p.value)
               ) %>%
        mutate(qvalue = p.adjust(pvalue, method = "BH"))
    ## perform wilcox test on each transcript
    return(out)
}

C08leaf_polyA_changed <- compare_polyA_between_samples(polyA_profile,
                                    "C08leaf_0h",
                                    "C08leaf_1h") %>%
    filter(qvalue <= 0.05) %>%
    mutate(deltaLength = sample2_median - sample1_median)

C08leaf_polyA_changed <- C08leaf_polyA_changed %>%
    select(transID, sample1_median, sample2_median,
           deltaLength, pvalue, qvalue) %>%
    rename(C08leaf_0h_polyA = sample1_median,
           C08leaf_1h_polyA = sample2_median)

C08root_polyA_changed <- compare_polyA_between_samples(polyA_profile,
                                                       "C08root_0h",
                                                       "C08root_1h") %>%
    filter(qvalue <= 0.05) %>%
    mutate(deltaLength = sample2_median - sample1_median)

C08root_polyA_changed <- C08root_polyA_changed %>%
    select(transID, sample1_median, sample2_median,
           deltaLength, pvalue, qvalue) %>%
    rename(C08root_0h_polyA = sample1_median,
           C08root_1h_polyA = sample2_median)

write.table(C08leaf_polyA_changed,
            file = "alternative_polyA_length.C08leaf.tsv",
            quote = F, sep = "\t", row.names = FALSE,
            col.names = TRUE)

write.table(C08root_polyA_changed,
            file = "alternative_polyA_length.C08root.tsv",
            quote = F, sep = "\t", row.names = FALSE,
            col.names = TRUE)
