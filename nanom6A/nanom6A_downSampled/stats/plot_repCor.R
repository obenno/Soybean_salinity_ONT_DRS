#! /usr/bin/env Rscript

library(tidyverse)
library(patchwork)
library(ggpubr)
library(foreach)

Leaf_0h_data <- read_tsv("Leaf_0h.ratioCommon.tsv") %>%
    filter(total_n.r1 >= 20,
           total_n.r2 >= 20,
           total_n.r3 >= 20) %>%
    select(geneID, chr, pos, motif, `ratio.r1`, `ratio.r2`, `ratio.r3`)

Leaf_0h_R1R2_cor <- cor(Leaf_0h_data$ratio.r1, Leaf_0h_data$ratio.r2)
Leaf_0h_R1R3_cor <- cor(Leaf_0h_data$ratio.r1, Leaf_0h_data$ratio.r3)
Leaf_0h_R2R3_cor <- cor(Leaf_0h_data$ratio.r2, Leaf_0h_data$ratio.r3)

Leaf_1h_data <- read_tsv("Leaf_1h.ratioCommon.tsv") %>%
    filter(total_n.r1 >= 20,
           total_n.r2 >= 20,
           total_n.r3 >= 20) %>%
    select(geneID, chr, pos, motif, `ratio.r1`, `ratio.r2`, `ratio.r3`)

Leaf_1h_R1R2_cor <- cor(Leaf_1h_data$ratio.r1, Leaf_1h_data$ratio.r2)
Leaf_1h_R1R3_cor <- cor(Leaf_1h_data$ratio.r1, Leaf_1h_data$ratio.r3)
Leaf_1h_R2R3_cor <- cor(Leaf_1h_data$ratio.r2, Leaf_1h_data$ratio.r3)

Root_0h_data <- read_tsv("Root_0h.ratioCommon.tsv") %>%
    filter(total_n.r1 >= 20,
           total_n.r2 >= 20,
           total_n.r3 >= 20) %>%
    select(geneID, chr, pos, motif, `ratio.r1`, `ratio.r2`, `ratio.r3`)

Root_0h_R1R2_cor <- cor(Root_0h_data$ratio.r1, Root_0h_data$ratio.r2)
Root_0h_R1R3_cor <- cor(Root_0h_data$ratio.r1, Root_0h_data$ratio.r3)
Root_0h_R2R3_cor <- cor(Root_0h_data$ratio.r2, Root_0h_data$ratio.r3)


Root_1h_data <- read_tsv("Root_1h.ratioCommon.tsv") %>%
    filter(total_n.r1 >= 20,
           total_n.r2 >= 20,
           total_n.r3 >= 20) %>%
    select(geneID, chr, pos, motif, `ratio.r1`, `ratio.r2`, `ratio.r3`)

Root_1h_R1R2_cor <- cor(Root_1h_data$ratio.r1, Root_1h_data$ratio.r2)
Root_1h_R1R3_cor <- cor(Root_1h_data$ratio.r1, Root_1h_data$ratio.r3)
Root_1h_R2R3_cor <- cor(Root_1h_data$ratio.r2, Root_1h_data$ratio.r3)


plotCor <- function(df, xcol, ycol, corNum, xlabel, ylabel){
    p <- ggplot(
        df,
        aes_string(x=xcol, y=ycol)) +
    geom_point(alpha = 0.6, color="#7570b3") +
    annotate("text", x = 0.05, y = 0.95,
             label = paste0("cor = ", signif(corNum,3)),
             hjust = 0) +
    annotate("text", x = 0.05, y = 0.85,
             label = paste0("n = ", nrow(df)),
             hjust = 0) +
    theme_bw() + xlab(xlabel) + ylab(ylabel) +
    theme(panel.grid = element_blank())
}

p1 <- plotCor(Leaf_0h_data, 'ratio.r1', 'ratio.r2', Leaf_0h_R1R2_cor, 'Leaf_0h_R1', 'Leaf_0h_R2')
p2 <- plotCor(Leaf_0h_data, 'ratio.r1', 'ratio.r3', Leaf_0h_R1R3_cor, 'Leaf_0h_R1', 'Leaf_0h_R3')
p3 <- plotCor(Leaf_0h_data, 'ratio.r2', 'ratio.r3', Leaf_0h_R2R3_cor, 'Leaf_0h_R2', 'Leaf_0h_R3')

p4 <- plotCor(Leaf_1h_data, 'ratio.r1', 'ratio.r2', Leaf_1h_R1R2_cor, 'Leaf_1h_R1', 'Leaf_1h_R2')
p5 <- plotCor(Leaf_1h_data, 'ratio.r1', 'ratio.r3', Leaf_1h_R1R3_cor, 'Leaf_1h_R1', 'Leaf_1h_R3')
p6 <- plotCor(Leaf_1h_data, 'ratio.r2', 'ratio.r3', Leaf_1h_R2R3_cor, 'Leaf_1h_R2', 'Leaf_1h_R3')

p7 <- plotCor(Root_0h_data, 'ratio.r1', 'ratio.r2', Root_0h_R1R2_cor, 'Root_0h_R1', 'Root_0h_R2')
p8 <- plotCor(Root_0h_data, 'ratio.r1', 'ratio.r3', Root_0h_R1R3_cor, 'Root_0h_R1', 'Root_0h_R3')
p9 <- plotCor(Root_0h_data, 'ratio.r2', 'ratio.r3', Root_0h_R2R3_cor, 'Root_0h_R2', 'Root_0h_R3')

p10 <- plotCor(Root_1h_data, 'ratio.r1', 'ratio.r2', Root_1h_R1R2_cor, 'Root_1h_R1', 'Root_1h_R2')
p11 <- plotCor(Root_1h_data, 'ratio.r1', 'ratio.r3', Root_1h_R1R3_cor, 'Root_1h_R1', 'Root_1h_R3')
p12 <- plotCor(Root_1h_data, 'ratio.r2', 'ratio.r3', Root_1h_R2R3_cor, 'Root_1h_R2', 'Root_1h_R3')

patchwork <- p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12 +  plot_layout(ncol = 3)

ggsave(patchwork, filename = "m6A_pairwise_cor.png", width = 10, height = 12)

Leaf_0h_transformed <- Leaf_0h_data %>%
    pivot_longer(
        c(`ratio.r1`, `ratio.r2`, `ratio.r3`),
        names_to = 'replicate',
        values_to = 'ratio'
    ) %>%
    mutate(sample="Leaf_0h")

Leaf_1h_transformed <- Leaf_1h_data %>%
    pivot_longer(
        c(`ratio.r1`, `ratio.r2`, `ratio.r3`),
        names_to = 'replicate',
        values_to = 'ratio'
    ) %>%
    mutate(sample="Leaf_1h")

Root_0h_transformed <- Root_0h_data %>%
    pivot_longer(
        c(`ratio.r1`, `ratio.r2`, `ratio.r3`),
        names_to = 'replicate',
        values_to = 'ratio'
    ) %>%
    mutate(sample="Root_0h")

Root_1h_transformed <- Root_1h_data %>%
    pivot_longer(
        c(`ratio.r1`, `ratio.r2`, `ratio.r3`),
        names_to = 'replicate',
        values_to = 'ratio'
    ) %>%
    mutate(sample="Root_1h")

## write m6A data to tables for supplementary,
## each tissue separately, not combined table
write_tsv(Leaf_0h_data, file = "m6A_ratio.Leaf_0h_combined.tsv")
write_tsv(Leaf_1h_data, file = "m6A_ratio.Leaf_1h_combined.tsv")
write_tsv(Root_0h_data, file = "m6A_ratio.Root_0h_combined.tsv")
write_tsv(Root_1h_data, file = "m6A_ratio.Root_1h_combined.tsv")

fullData <- bind_rows(Leaf_0h_transformed,
                      Leaf_1h_transformed,
                      Root_0h_transformed,
                      Root_1h_transformed) %>%
    mutate(replicate = str_replace(replicate, "ratio.r", "rep"))

##test_Res <- wilcox.test(Nodule_avg$ratio, Root_avg$ratio)$p.value %>%
##                                                        signif(3)

p <- ggplot(fullData,
             aes(x=sample, y=ratio, fill = replicate)) +
    geom_boxplot(color = "black", alpha = 0.8) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_brewer(palette="Dark2")## +
    ##annotate("text", x = 1.5, y = 0.8,
    ##         label = paste0("p = NoP"))

ggsave(p, filename = "m6A_ratio.boxplot.pdf", width = 5, height = 5)

fullData_ratioMean <- fullData %>%
    group_by(geneID, chr, pos, motif, sample) %>%
    summarise(ratioMean=mean(ratio))

saveRDS(fullData_ratioMean, "fullData_ratioMean.rds")

p_ratioMean <- ggplot(fullData_ratioMean,
                      aes(x=sample, y=ratioMean,
                          fill = sample))+
    geom_boxplot(color = "black", alpha = 0.8) +
    ylab("m6A ratio") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_brewer(palette="Dark2")

my_comparisons <- list(c("Leaf_0h", "Leaf_1h"),
                       c("Root_0h", "Root_1h"),
                       c("Leaf_0h", "Root_0h"),
                       c("Leaf_1h", "Root_1h"))

p_ratioMean <- p_ratioMean +
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       tip.length = 0)

ggsave(plot=p_ratioMean,
       filename = "m6A_ratio.boxplot.ratioMean.pdf",
       width = 5, height = 5)

leaf_up_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Leaf_1hvs0h.gene_up.out") %>% pull(geneID)

leaf_down_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Leaf_1hvs0h.gene_down.out") %>% pull(geneID)

root_up_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Root_1hvs0h.gene_up.out") %>% pull(geneID)

root_down_genes <- read_tsv("~/salinity_suppl_analysis/analysis/expressions/Root_1hvs0h.gene_down.out") %>% pull(geneID)

m6Aratio_catDEG <- fullData_ratioMean %>%
    mutate(
        DEG=case_when(
            geneID %in% leaf_up_genes & sample %in% c("Leaf_0h", "Leaf_1h") ~ "Leaf_up_DEG",
            geneID %in% leaf_down_genes & sample %in% c("Leaf_0h", "Leaf_1h") ~ "Leaf_down_DEG",
            geneID %in% root_up_genes & sample %in% c("Root_0h", "Root_1h") ~ "Root_up_DEG",
            geneID %in% root_down_genes & sample %in% c("Root_0h", "Root_1h") ~ "Root_down_DEG",
            str_detect(sample, "Leaf") ~ "Leaf_non_DEG",
            str_detect(sample, "Root") ~ "Root_non_DEG",
            TRUE ~ "noCategory"
        )
    ) %>%
    mutate(DEG = factor(DEG, levels = c("Leaf_up_DEG",
                                        "Leaf_non_DEG",
                                        "Leaf_down_DEG",
                                        "Root_up_DEG",
                                        "Root_non_DEG",
                                        "Root_down_DEG")))

## Get mann-whitney test p
##m6Aratio_catDEG_p <- m6Aratio_catDEG %>%
##    group_by(DEG) %>%
##    rstatix::wilcox_test(ratioMean ~ sample)

p_m6A_DEG <- ggplot(m6Aratio_catDEG,
                    aes(x=sample, y=ratioMean, fill=DEG)) +
    geom_boxplot(color = "black", alpha = 0.8) +
    ylab("m6A ratio") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))+
    scale_fill_brewer(palette="Dark2")

## Compute Kruskal-Wallis test p value for each time point
samples <- c("Leaf_0h", "Leaf_1h", "Root_0h", "Root_1h")
kw_pvalue <- foreach(i=1:4, .combine = c) %do% {
    t <- m6Aratio_catDEG %>%
        filter(sample==samples[i])
    kruskal.test(ratioMean ~ DEG, data = t)$p.value
}


kw_pvalue <- tibble(pvalue=kw_pvalue) %>%
    mutate(
        pvalue_label=case_when(
            pvalue<0.001 ~ scales::label_scientific(digits = 3)(pvalue),
            TRUE ~ scales::label_number(0.001)(pvalue)
        )
    ) %>% pull(pvalue_label)


##p_m6A_DEG <- p_m6A_DEG +
##    geom_text(data = m6Aratio_catDEG_p,
##              aes(x=DEG, y=1.05, label=p, fill =NA))
    ##stat_pvalue_manual(stat.test,  label = "p", tip.length = 0)

p_m6A_DEG <- p_m6A_DEG +
    ggplot2::annotate("text", x=1:4, y=1.05,
                      label=kw_pvalue)

ggsave(p_m6A_DEG,
       filename = "m6A_ratio_in_DEG.boxplot.pdf",
       width = 6,
       height = 5)


## create m6A vs DEG with pair-wise mann-whitney test
m6A_mw_stat <- m6Aratio_catDEG %>%
    group_by(sample) %>%
    rstatix::wilcox_test(ratioMean~DEG, p.adjust.method = "none")

m6A_mw_stat <- m6A_mw_stat %>%
    filter(str_detect(group1, "non") | str_detect(group2,"non")) %>%
    mutate(
        y_position = case_when(
            str_detect(group2, "non") ~ 1.05,
            TRUE ~ 1.15
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
##x_min=c(c(1,2,3,4)-0.25, c(1,2,3,4)+0.05),


p_m6A_DEG <- ggplot(m6Aratio_catDEG,
                    aes(x=sample, y=ratioMean, fill=DEG)) +
    geom_boxplot(color = "black", alpha = 0.8) +
    ylab("m6A ratio") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))+
    scale_fill_brewer(palette="Dark2")

p_m6A_DEG <- p_m6A_DEG +
    geom_text(data = m6A_mw_stat,
              aes(x=sample, y=y_position, group=`group`,
                  label = `p.adj.signif`),
              inherit.aes = FALSE,
              position = position_dodge(width = 0.6))+
    geom_segment(data = m6A_mw_stat,
                 aes(x=x_min,
                     xend=x_max,
                     y=y_position-0.03,
                     yend=y_position-0.03),
                 color="black",
                 inherit.aes = FALSE) +
    theme(legend.position = "bottom")

ggsave(p_m6A_DEG,
       filename = "m6A_ratio_in_DEG.boxplot.mwTest.pdf",
       width = 6,
       height = 5)


## create DEG/non-DEG m6A ratio sample 0h vs 1h with pair-wise mann-whitney test
m6A_mw_stat <- m6Aratio_catDEG %>%
    group_by(DEG) %>%
    rstatix::wilcox_test(ratioMean~sample, p.adjust.method = "none")

m6A_mw_stat <- m6A_mw_stat %>%
    mutate(
        x_position = c(1,2,3,4,5,6),
        x_min=x_position-0.25,
        x_max=x_position+0.25
    )

m6A_DEG_expr_n <- m6Aratio_catDEG %>%
    group_by(DEG, sample) %>%
    dplyr::filter(!is.na(ratioMean)) %>%
    summarise(PAC_count = n()) %>%
    mutate(PAC_count = scales::label_comma()(PAC_count)) %>%
    mutate(
        y_position = case_when(
            str_detect(sample, "0h") ~ 1.05,
            TRUE ~ 1.1
        )
    )


p_m6A_DEG <- ggplot(m6Aratio_catDEG,
                    aes(x=DEG, y=ratioMean, fill=sample)) +
    geom_boxplot(color = "black", alpha = 0.8) +
    ylab("m6A ratio") + xlab("DEGs group")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom")+
    scale_fill_brewer(palette="Dark2")

p_m6A_DEG <- p_m6A_DEG +
    geom_text(data = m6A_DEG_expr_n,
              aes(DEG, y=y_position, group=sample,
                  label = `PAC_count`),
              inherit.aes = FALSE,
              position = position_dodge(width = 0.9)) +
    geom_text(data = m6A_mw_stat,
              aes(x=DEG, y=1.2,
                  label = `p.adj.signif`),
              inherit.aes = FALSE) +
    geom_segment(data = m6A_mw_stat,
                 aes(x=x_min,
                     xend=x_max,
                     y=1.15,
                     yend=1.15),
                 color="black",
                 inherit.aes = FALSE)

ggsave(p_m6A_DEG,
       filename = "m6A_ratio_in_DEG.boxplot.mwTest.byDEG.pdf",
       width = 7,
       height = 5)

## Generate scatter plot for the m6A sites detected in both 0h and 1h
m6A_sites_avg <- fullData %>%
    group_by(geneID, chr, pos, motif, sample) %>%
    summarise(ratio = mean(ratio)) %>%
    pivot_wider(c(geneID, chr, pos, motif), names_from = "sample", values_from = "ratio")


plotScatter <- function (df, colx, coly) {
    d <- df %>%
        select(`geneID`, `chr`, `pos`, `motif`, all_of(c(colx, coly))) %>%
        na.omit %>%
        mutate(change = case_when(.data[[coly]]- .data[[colx]] >= 0.1 ~ "up",
                                  .data[[coly]]- .data[[colx]] <=- 0.1 ~ "down",
                                  TRUE ~ "unchanged"))
    up_number <- d %>% filter(change == "up") %>% nrow()
    down_number <- d %>% filter(change == "down") %>% nrow()
    unchanged_number <- d %>% filter(change == "unchanged") %>% nrow()

    p <- ggplot(d, aes_string(x=colx, y=coly, color="change")) +
        geom_point(alpha = 0.6) +
        annotate(
            "text", x = 0.05, y= 0.9,
            label = paste0("up: ", up_number, "\n",
                           "unchanged: ", unchanged_number, "\n",
                           "down: ", down_number),
            hjust = 0
        ) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "none") +
        scale_color_manual(values = c("#0571b0", "grey", "#ca0020"))
}

p1 <- plotScatter(m6A_sites_avg, "Leaf_0h", "Leaf_1h")
p2 <- plotScatter(m6A_sites_avg, "Root_0h", "Root_1h")
p3 <- plotScatter(m6A_sites_avg, "Leaf_0h", "Root_0h")
p4 <- plotScatter(m6A_sites_avg, "Leaf_1h", "Root_1h")

patchwork1 <- p1+p2
ggsave(patchwork1, filename = "m6A_ratio.scatter.pdf", width = 10, height = 5)

patchwork2 <- p3+p4
ggsave(patchwork2, filename = "m6A_ratio.leaf_vs_root.scatter.pdf", width = 10, height = 5)

## Plot delta_m6A vs geneFC
m6Aratio_catDEG_deltaRatio <- m6Aratio_catDEG %>%
    ##select(-DEG) %>%
    pivot_wider(c(geneID, chr, pos, motif, DEG),
                names_from = "sample",
                values_from = "ratioMean") %>%
    mutate(deltaLeaf=Leaf_1h-Leaf_0h, deltaRoot=Root_1h-Root_0h) %>%
    filter(!str_detect(DEG, "non")) %>%
    select(geneID, chr, pos, motif, deltaLeaf, deltaRoot)

LeafFC <- read_tsv("/home/ubuntu/salinity_suppl_analysis/analysis/expressions/Leaf_1hvs0h.DEGresult.out") %>%
    select(geneID, log2FoldChange)

RootFC <- read_tsv("/home/ubuntu/salinity_suppl_analysis/analysis/expressions/Root_1hvs0h.DEGresult.out") %>%
    select(geneID, log2FoldChange)

m6Aratio_catDEG_deltaRatio <- m6Aratio_catDEG_deltaRatio %>%
    left_join(LeafFC) %>%
    dplyr::rename("Leaf_logFC"=log2FoldChange) %>%
    left_join(RootFC) %>%
    dplyr::rename("Root_logFC"=log2FoldChange)

m6Aratio_catDEG_deltaRatio_leaf <- m6Aratio_catDEG_deltaRatio %>%
    select(-c(deltaRoot, Root_logFC)) %>%
    mutate(
        category=case_when(
            deltaLeaf >= 0.1 & Leaf_logFC >= 1 ~ "Leaf_positive_cor",
            deltaLeaf >= 0.1 & Leaf_logFC <= -1 ~ "Leaf_negative_cor",
            deltaLeaf <= -0.1 & Leaf_logFC >= 1 ~ "Leaf_negative_cor",
            deltaLeaf <= -0.1 & Leaf_logFC <= -1 ~ "Leaf_positive_cor",
            TRUE ~ "none"
            )
    )

m6Aratio_catDEG_deltaRatio_root <- m6Aratio_catDEG_deltaRatio %>%
    select(-c(deltaLeaf, Leaf_logFC)) %>%
    mutate(
        category=case_when(
            deltaRoot >= 0.1 & Root_logFC >= 1 ~ "Root_positive_cor",
            deltaRoot >= 0.1 & Root_logFC <= -1 ~ "Root_negative_cor",
            deltaRoot <= -0.1 & Root_logFC >= 1 ~ "Root_negative_cor",
            deltaRoot <= -0.1 & Root_logFC <= -1 ~ "Root_positive_cor",
            TRUE ~ "none"
            )
    )


p1 <- ggplot(m6Aratio_catDEG_deltaRatio_leaf,
             aes(deltaLeaf, Leaf_logFC, color=category))+
    geom_point(alpha=0.6)+
    xlim(-0.5,0.5)+
    ylim(-4,4)+
    scale_color_manual(values = c("Leaf_negative_cor" = "#0571b0",
                                  "none"="grey",
                                  "Leaf_positive_cor" = "#ca0020"))+
    theme_bw()+theme(panel.grid = element_blank(),
                     legend.position = "None")
p2 <- ggplot(m6Aratio_catDEG_deltaRatio_root,
             aes(deltaRoot, Root_logFC, color=category))+
    geom_point(alpha=0.6)+
    xlim(-0.5,0.5)+
    ylim(-4,4)+
    scale_color_manual(values = c("Root_negative_cor" = "#0571b0",
                                  "none"="grey",
                                  "Root_positive_cor" = "#ca0020"))+
    theme_bw()+theme(panel.grid = element_blank(),
                     legend.position = "None")
p <- p1+p2
ggsave("deltaRatio_vs_logFC.pdf", width = 6, height = 3)

m6Aratio_catDEG_deltaRatio_leaf %>%
    filter(Leaf_logFC>=1|Leaf_logFC<=-1,
           deltaLeaf>=0.1|deltaLeaf<=-0.1) %>%
    write_tsv(file = "deltaRatio_vs_logFC.leaf.tsv")

m6Aratio_catDEG_deltaRatio_root %>%
    filter(Root_logFC>=1|Root_logFC<=-1,
           deltaRoot>=0.1|deltaRoot<=-0.1) %>%
    write_tsv(file = "deltaRatio_vs_logFC.root.tsv")
