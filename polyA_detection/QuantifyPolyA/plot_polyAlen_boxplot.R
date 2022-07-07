#! /usr/bin/env Rscript

library(tidyverse)
library(foreach)
library(ggpubr)

files <- dir(pattern = ".polyALen.tsv")

d <- foreach(i=1:length(files), .combine = bind_rows) %do% {
    sampleName <- str_replace(files[i], ".polyALen.tsv", "")
    read_tsv(files[i], col_names=c("chr", "start", "end", "polyA", "score", "strand", "len")) %>%
        mutate(sample = sampleName)
}

d <- d %>%
    mutate(tissue=str_replace(sample, "_R[123]$", ""))

p <- ggplot(d, aes(x=sample, y=len, fill=tissue)) +
    geom_violin(alpha=0.6)+
    geom_boxplot(alpha=0.6, width=0.2, outlier.shape = NA)+
    theme_bw()+
    xlab("Library") + ylab("polyA Length (nt)") +
    scale_fill_brewer(palette = "Dark2")+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))
ggplot2::ggsave(p, filename = "polyA_len_boxplot.pdf",
                width = 5, height = 5)

polyAlen_mean <- d %>%
    group_by(chr, start, end, polyA, score, strand, tissue) %>%
    summarise(len=mean(len)) # don't use na.rm=TRUE since 3 replicates are required

p <- ggplot(polyAlen_mean,
            aes(x=tissue, y=len, fill=tissue)) +
    geom_violin(alpha = 0.6)+
    geom_boxplot(alpha=0.6, width=0.2, outlier.shape = NA)+
    theme_bw()+
    xlab("Sample") + ylab("polyA Length (nt)") +
    scale_fill_brewer(palette = "Dark2")+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))

my_comparisons <- list(c("Leaf_0h", "Leaf_1h"),
                       c("Root_0h", "Root_1h"),
                       c("Leaf_0h", "Root_0h"),
                       c("Leaf_1h", "Root_1h"))

p <- p + stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       tip.length = 0)

ggplot2::ggsave(p, filename = "polyA_len_replicateMean.boxplot.pdf",
                 width = 5, height = 5)

polyA_info <- d %>% dplyr::select(-tissue) %>%
    pivot_wider(names_from = `sample`, values_from = `len`) %>%
    mutate(geneID=str_extract(polyA, "XLOC_[0-9]{6}"),
           totalCount = str_extract(polyA, "[0-9]+$"))

## select top two polyA sites per gene loci by totalCount
## add type infor of proximal or distal
polyA_proximal_distal <- polyA_info %>%
    dplyr::filter(!str_detect(polyA, "intergenic")) %>%
    group_by(geneID) %>%
    dplyr::filter(n()>=2) %>%
    arrange(desc(totalCount), .by_group = TRUE) %>%
    slice_head(n=2) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(
        type=case_when(
            strand == "+" & row_number() ==1 ~ "proximal",
            strand == "+" & row_number() ==2 ~ "distal",
            strand == "-" & row_number() ==1 ~ "distal",
            strand == "-" & row_number() ==2 ~ "proximal",
        )
    ) %>%
    ungroup()

write_tsv(polyA_proximal_distal, file = "polyA_proximal_distal.info.tsv")

polyA_proximal_distal <- polyA_proximal_distal %>%
    pivot_longer(c(starts_with(c("Leaf", "Root"))),
                 names_to = "sample", values_to = "len")

p <- ggplot(polyA_proximal_distal, aes(x=sample, y=len, fill=type))+
    geom_violin(alpha=0.6)+
    geom_boxplot(alpha=0.8, width=0.2, position = position_dodge(width = 0.9),
                 na.rm=TRUE, show.legend = FALSE, outlier.shape = NA)+
    theme_bw()+
    xlab("Library") + ylab("polyA Length (nt)") +
    scale_fill_brewer(palette = "Dark2")+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))
ggplot2::ggsave(p, filename = "polyA_len_proximal_distal_boxplot.pdf")

## investigate polyA length vs gene expression

## calculate mean polyA length of proximal and distal
polyA_proximal_distal <- polyA_proximal_distal %>%
    rename("library" = sample) %>%
    mutate(sample = str_replace(library, "_R[123]$", "")) %>%
    group_by(chr, start, end, polyA,
             score, strand, geneID, totalCount,
             type, sample) %>%
    summarise(polyA_len=mean(len, na.rm = TRUE))

## get single PAC gene infor
polyA_single <- polyA_info %>%
    dplyr::filter(!str_detect(polyA, "intergenic")) %>%
    group_by(geneID) %>%
    dplyr::filter(n()==1) %>%
    mutate(type = "single") %>%
    ungroup()

polyA_single <- polyA_single %>%
    pivot_longer(c(starts_with(c("Leaf", "Root"))),
                 names_to = "library", values_to = "len") %>%
    mutate(sample = str_replace(library, "_R[123]$", "")) %>%
    group_by(chr, start, end, polyA,
             score, strand, geneID, totalCount,
             type, sample) %>%
    summarise(polyA_len=mean(len, na.rm = TRUE))


polyA_dataset <- bind_rows(polyA_single,
                           polyA_proximal_distal) %>%
    ungroup()

geneTPM <- read_tsv("/home/ubuntu/salinity_suppl_analysis/analysis/expressions/gene.tpm.tsv") %>%
    dplyr::select(-ref) %>%
    pivot_longer(-geneID,
                 names_to = "library",
                 values_to = "TPM") %>%
    mutate(sample=str_replace(library, "_R[123]$", "")) %>%
    group_by(geneID, sample) %>%
    summarise(TPM=mean(TPM, na.rm = TRUE))

polyA_vs_geneExpr <- polyA_dataset %>%
    left_join(geneTPM, by = c("geneID", "sample"))

## calculate correlation
cor_tbl <- polyA_vs_geneExpr %>%
    group_by(type, sample) %>%
    summarise(correlation = cor(polyA_len, TPM,
                                method = "spearman",
                                use = "complete.obs")) %>%
    mutate(cor_label = paste0("Cor = ",
                                scales::label_number(0.001)(correlation)))

p <- ggplot(polyA_vs_geneExpr,
            aes(x=polyA_len, y=log2(TPM+1), color=type))+
    geom_point(alpha=0.4) +
    scale_color_brewer(palette = "Dark2")+
    geom_text(data=cor_tbl,
              aes(x=Inf, y=Inf,
                  hjust = 1.1, vjust = 2,
                  label = cor_label),
              show.legend = FALSE,
              inherit.aes = FALSE) +
    facet_grid(type~sample)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank())

ggsave("polyAlen_vs_geneExpr.scatter.png",
       width = 7, height = 6)

## get DEG list
leaf_up_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Leaf_1hvs0h.gene_up.out") %>% pull(geneID)

leaf_down_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Leaf_1hvs0h.gene_down.out") %>% pull(geneID)

root_up_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Root_1hvs0h.gene_up.out") %>% pull(geneID)

root_down_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Root_1hvs0h.gene_down.out") %>% pull(geneID)

## check polyA tail length vs TPM of DEGs
polyA_vs_DEG_expr <- polyA_vs_geneExpr %>%
    mutate(
        DEG = case_when(
            geneID %in% leaf_up_genes & str_detect(sample, "Leaf")  ~ "Leaf_up_DEG",
            geneID %in% leaf_down_genes & str_detect(sample, "Leaf") ~ "Leaf_down_DEG",
            geneID %in% root_up_genes & str_detect(sample, "Root") ~ "Root_up_DEG",
            geneID %in% root_down_genes & str_detect(sample, "Root") ~ "Root_down_DEG",
            str_detect(sample, "Leaf") ~ "Leaf_non_DEG",
            str_detect(sample, "Root") ~ "Root_non_DEG"
        )
    ) ##%>%
    ##dplyr::filter(DEG!="nonDEG")

polyA_vs_DEG_expr <- polyA_vs_DEG_expr %>%
    mutate(DEG=factor(DEG, levels = c("Leaf_up_DEG",
                                      "Leaf_non_DEG",
                                      "Leaf_down_DEG",
                                      "Root_up_DEG",
                                      "Root_non_DEG",
                                      "Root_down_DEG")))

polyA_vs_DEG_expr_n <- polyA_vs_DEG_expr %>%
    group_by(sample, DEG) %>%
    dplyr::filter(!is.na(polyA_len)) %>%
    summarise(PAC_count = n()) %>%
    mutate(PAC_count = scales::label_comma()(PAC_count)) %>%
    mutate(
        y_position = case_when(
            str_detect(DEG, "non") ~ 460,
            TRUE ~ 440
        )
    )


polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr %>%
    group_by(sample) %>%
    rstatix::wilcox_test(polyA_len~DEG, p.adjust.method = "none")

polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr_p %>%
    filter(str_detect(group1, "non") | str_detect(group2,"non")) %>%
    mutate(
        y_position = case_when(
            str_detect(group2, "non") ~ 400,
            TRUE ~ 380
        ),
        ##group = paste0(group1, "_vs_", group2)
        group = row_number()
    ) %>%
    mutate(
        x_min=rep(c(1,2,3,4), each = 2),
        x_min=case_when(
            row_number() %% 2 ==1 ~ x_min-0.25,
            row_number() %% 2 ==0 ~ x_min+0.05
        ),
        x_max=x_min+0.2
    )


p <- ggplot(polyA_vs_DEG_expr,
            aes(x=sample, y=polyA_len, fill=DEG))+
    geom_boxplot(color = "black", alpha=0.8) +
    scale_fill_brewer(palette = "Dark2")+
    xlab("DEGs group") + ylab("polyA length (nt) of PACs") +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank())

p <- p + geom_text(data = polyA_vs_DEG_expr_n,
                   aes(y=y_position, label = PAC_count, group = DEG),
                   position = position_dodge(width = .6))

p <- p + geom_text(data = polyA_vs_DEG_expr_p,
                   aes(x= sample, y= y_position, label = `p.adj.signif`,
                       group = group),
                   inherit.aes = FALSE,
                   position = position_dodge(width = 0.6))

p <- p + geom_segment(data = polyA_vs_DEG_expr_p,
                      aes(x=x_min,
                          xend=x_max,
                          y=y_position-5,
                          yend=y_position-5),
                      color="black",
                      inherit.aes = FALSE)

ggsave("polyAlen_vs_DEG_expr.boxplot.pdf",
       width = 7, height = 6)

## Plot polyA length vs DEG
## (but multiple PACs on the same gene will be averaged)
polyA_vs_DEG_expr <- polyA_vs_DEG_expr %>%
    group_by(sample, DEG, geneID) %>%
    summarise(polyA_len = mean(polyA_len, na.rm = TRUE))

polyA_vs_DEG_expr_n <- polyA_vs_DEG_expr %>%
    group_by(sample, DEG) %>%
    dplyr::filter(!is.na(polyA_len)) %>%
    summarise(PAC_count = n()) %>%
    mutate(PAC_count = scales::label_comma()(PAC_count)) %>%
    mutate(
        y_position = case_when(
            str_detect(DEG, "non") ~ 460,
            TRUE ~ 440
        )
    )


polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr %>%
    group_by(sample) %>%
    rstatix::wilcox_test(polyA_len ~ DEG)

polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr_p %>%
    filter(str_detect(group1, "non") | str_detect(group2,"non")) %>%
    mutate(
        y_position = case_when(
            str_detect(group2, "non") ~ 400,
            TRUE ~ 380
        ),
        ##group = paste0(group1, "_vs_", group2)
        group = row_number()
    ) %>%
    mutate(
        x_min=rep(c(1,2,3,4), each = 2),
        x_min=case_when(
            row_number() %% 2 ==1 ~ x_min-0.25,
            row_number() %% 2 ==0 ~ x_min+0.05
        ),
        x_max=x_min+0.2
    )

p <- ggplot(polyA_vs_DEG_expr,
            aes(x=sample, y=polyA_len, fill=DEG))+
    geom_boxplot(color = "black", alpha=0.8) +
    scale_fill_brewer(palette = "Dark2")+
    xlab("DEGs group") + ylab("polyA length (nt) of PACs") +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank())

p <- p + geom_text(data = polyA_vs_DEG_expr_n,
                   aes(y=y_position, label = PAC_count, group = DEG),
                   position = position_dodge(width = .6))

p <- p + geom_text(data = polyA_vs_DEG_expr_p,
                   aes(x= sample, y= y_position, label = `p.adj.signif`,
                       group = group),
                   inherit.aes = FALSE,
                   position = position_dodge(width = 0.6))

p <- p + geom_segment(data = polyA_vs_DEG_expr_p,
                      aes(x=x_min,
                          xend=x_max,
                          y=y_position-5,
                          yend=y_position-5),
                      color="black",
                      inherit.aes = FALSE)

ggsave("geneMean_polyAlen_vs_DEG_expr.boxplot.pdf",
       width = 7, height = 6)


polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr %>%
    group_by(DEG) %>%
    rstatix::wilcox_test(polyA_len ~ sample)

polyA_vs_DEG_expr_p <- polyA_vs_DEG_expr_p %>%
    mutate(
        x_position = c(1,2,3,4,5,6),
        x_min=x_position-0.25,
        x_max=x_position+0.25
    )

polyA_vs_DEG_expr_n <- polyA_vs_DEG_expr %>%
    group_by(DEG, sample) %>%
    dplyr::filter(!is.na(polyA_len)) %>%
    summarise(PAC_count = n()) %>%
    mutate(PAC_count = scales::label_comma()(PAC_count)) %>%
    mutate(
        y_position = case_when(
            str_detect(sample, "0h") ~ 410,
            TRUE ~ 430
        )
    )

p <- ggplot(polyA_vs_DEG_expr,
            aes(x=DEG, y=polyA_len, fill=sample)) +
    geom_boxplot(color = "black", alpha=0.8) +
    scale_fill_brewer(palette = "Dark2")+
    xlab("DEGs group") + ylab("polyA length (nt) of PACs") +
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank())

p <- p + geom_text(data=polyA_vs_DEG_expr_n,
                   aes(x=DEG, y=y_position, group = sample,
                       label = PAC_count),
                   inherit.aes = FALSE,
                   position = position_dodge(0.9))

p <- p + geom_text(data=polyA_vs_DEG_expr_p,
                   aes(x=`x_position`, y=460, label = `p.adj.signif`),
                   inherit.aes = FALSE)

p <- p + geom_segment(data=polyA_vs_DEG_expr_p,
                      aes(x=`x_min`, y=450, xend=`x_max`, yend=450),
                      inherit.aes = FALSE)

ggsave("geneMean_polyAlen_vs_DEG_expr.byDEG.boxplot.pdf",
       width = 7, height = 6)