#! /usr/bin/env Rscript


set.seed(123)

library(tidyverse)

event_expr <- read_tsv("AS_aPAS_event_expr.tsv",
                       col_names = c('geneID', 'type',
                                     'inclusionExon', 'exclusionExon',
                                     'inclusionExpr', 'exclusionExpr'))

samples = c("C08leaf_0h", "C08leaf_1h", "C08root_0h", "C08root_1h")

calc_PSI <- function(inclusionExpr_string,
                     exclusionExpr_string){
    inclusion_expr <- unlist(str_split(inclusionExpr_string, fixed("|")))
    exclusion_expr <- unlist(str_split(exclusionExpr_string, fixed("|")))
    inclusion_expr <- map_dbl(inclusion_expr, as.numeric)
    exclusion_expr <- map_dbl(exclusion_expr, as.numeric)
    total_expr <- inclusion_expr + exclusion_expr
    psi <- inclusion_expr/total_expr
    psi_string <- paste(psi, collapse = "|")
    return(psi_string)
}

calc_deltaPSI <- function(inclusionExpr_string,
                     exclusionExpr_string){
    inclusion_expr <- unlist(str_split(inclusionExpr_string, fixed("|")))
    exclusion_expr <- unlist(str_split(exclusionExpr_string, fixed("|")))
    inclusion_expr <- map_dbl(inclusion_expr, as.numeric)
    exclusion_expr <- map_dbl(exclusion_expr, as.numeric)
    total_expr <- inclusion_expr + exclusion_expr
    psi <- inclusion_expr/total_expr
    ## Check ../alternative_splicing/format_AS_events_expr.awk
    ## for sample order
    if(!is.na(psi[2]) & !is.na(psi[1])){
        deltaPSI_1 <- psi[2] - psi[1]
    }else{
        deltaPSI_1 <- NA
    }
    if(!is.na(psi[4]) & !is.na(psi[3])){
        deltaPSI_2 <- psi[4] - psi[3]
    }else{
        deltaPSI_2 <- NA
    }
    return(c(deltaPSI_1, deltaPSI_2))
}

calc_fisherP <- function(inclusionExpr_string,
                         exclusionExpr_string){
    inclusion_expr <- unlist(str_split(inclusionExpr_string, fixed("|")))
    exclusion_expr <- unlist(str_split(exclusionExpr_string, fixed("|")))
    inclusion_expr <- map_dbl(inclusion_expr, as.numeric)
    exclusion_expr <- map_dbl(exclusion_expr, as.numeric)
    x=c(inclusion_expr[1], inclusion_expr[2])
    y=c(exclusion_expr[1], exclusion_expr[2])
    d=data.frame(x,y)
    p1 = fisher.test(d)$p.value
    x=c(inclusion_expr[3], inclusion_expr[4])
    y=c(exclusion_expr[3], exclusion_expr[4])
    d=data.frame(x,y)
    p2 = fisher.test(d)$p.value
    return(c(p1,p2))
}

event_expr <- event_expr %>%
    mutate(PSI=map2_chr(.x=inclusionExpr,
                        .y=exclusionExpr,
                        .f=calc_PSI))

event_expr <- event_expr %>%
    mutate(deltaPSI=map2(.x=inclusionExpr,
                         .y=exclusionExpr,
                         .f=calc_deltaPSI))

eventNum <- nrow(event_expr)
event_expr <- event_expr %>%
    unnest(deltaPSI) %>%
    mutate(sample = rep(c("C08leaf_deltaPSI",
                          "C08root_deltaPSI"),
                        eventNum)) %>%
    pivot_wider(names_from = sample, values_from = deltaPSI)

event_expr <- event_expr %>%
    mutate(Pvalue = map2(.x=inclusionExpr,
                         .y=exclusionExpr,
                         .f=calc_fisherP))

event_expr <- event_expr %>%
    unnest(Pvalue) %>%
    mutate(sample = rep(c("C08leaf_pvalue",
                          "C08root_pvalue"),
                        eventNum)) %>%
    pivot_wider(names_from = sample, values_from = Pvalue)

## Conduct FDR correction
event_expr  <- event_expr %>%
    mutate(C08leaf_FDR = p.adjust(C08leaf_pvalue, method = "BH"),
           C08root_FDR = p.adjust(C08root_pvalue, method = "BH"))

## write a full table fisrt
write.table(event_expr,
            file = "AS_aPAS_event_fullData.tsv",
            quote = F, sep = "\t",
            row.names = F,
            col.names = T)


dAS_C08leaf <- event_expr %>%
    filter(!is.na(C08leaf_deltaPSI),
           abs(C08leaf_deltaPSI)>=0.1,
           C08leaf_pvalue<=0.05) %>%
           ##C08leaf_FDR<=0.05) %>%
    select(geneID, type, inclusionExon, exclusionExon,
           C08leaf_deltaPSI, C08leaf_pvalue) %>%
    rename(deltaPSI = C08leaf_deltaPSI,
           Pvalue = C08leaf_pvalue) %>%
    mutate(comparison = "C08leaf_1hvs0h")

dAS_C08root <- event_expr %>%
    filter(!is.na(C08root_deltaPSI),
           abs(C08root_deltaPSI)>=0.1,
           C08root_pvalue<=0.05) %>%
           ##C08root_FDR<=0.05) %>%
    select(geneID, type, inclusionExon, exclusionExon,
           C08root_deltaPSI, C08root_pvalue) %>%
    rename(deltaPSI = C08root_deltaPSI,
           Pvalue = C08root_pvalue) %>%
    mutate(comparison = "C08root_1hvs0h")

dAS <- rbind(dAS_C08leaf, dAS_C08root)

write.table(dAS,
            file = "dAS_data.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

## Extract differential PAS
dPAS <- dAS %>%
    filter(type == "aPAS")

write.table(dPAS,
            file = "dPAS.tsv",
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

dAS_count  <- dAS %>%
    group_by(type, comparison) %>%
    summarize(count = n())

dAS_eventCount <- dAS %>%
    select(geneID, inclusionExon, exclusionExon, type) %>%
    unique() %>%
    group_by(type) %>%
    summarise(count = n())

dAS_geneCount <- dAS %>%
    select(geneID, type) %>%
    unique() %>%
    group_by(type) %>%
    summarize(count = n())

dAS_geneList <- dAS %>%
    select(geneID) %>%
    unique()

write.table(dAS_geneList, file = "DASG.lst",
            quote = FALSE, sep = "\t",
            row.names=FALSE, col.names=FALSE)

write.table(dAS_eventCount, file = "dAS_eventCount.tsv",
            quote = FALSE, sep = "\t",
            row.names=FALSE, col.names=TRUE)

write.table(dAS_geneCount, file = "dAS_geneCount.tsv",
            quote = FALSE, sep = "\t",
            row.names=FALSE, col.names=TRUE)

p <- ggplot(dAS, aes(x=type, y=deltaPSI, fill = type))
p <- p + geom_boxplot(alpha = 0.6)
##p <- p + geom_boxplot()
p <- p + geom_hline(yintercept=0, linetype = 'dashed',
                    color = 'black', alpha = 0.8)
p <- p + geom_text(data = dAS_count,
                   aes(x=type, y=1.1,
                       label = count))

p <- p + facet_grid(~comparison)
p <- p + ylab(expression(paste(Delta, "PSI", sep="")))
p <- p + xlab("Type of events")
p <- p + scale_fill_brewer(palette = 'Set1')
p <- p + theme_bw(base_size=15)
p <- p + theme(legend.position = "none",
               panel.grid = element_blank(),
               panel.border = element_rect(colour = "black",
                                           fill=NA, size = 0.6),
               axis.ticks = element_line(size = 0.6),
               axis.text = element_text(color = "black"),
               axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
p <- p + theme(strip.background = element_blank())
##p <- p + geom_violin()
##p <- p + geom_jitter(alpha = 0.5)


ggsave(filename = 'PSI.pdf', width = 6, height = 4.5)


