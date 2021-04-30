#! /usr/bin/env Rscript

## This script will help to generate raw aPAS data
## for delta PSI calculation and differential test

library(tidyverse)


pas <- read_tsv("collapsed_gene_PAS_reads.strand.countsTable.tsv",
              col_names = c("geneID",
                            "chr", "start", "end", "strand",
                            "C08leaf_0h", "C08leaf_1h",
                            "C08root_0h", "C08root_1h"))
## Filter with max count >=5
pas <- pas %>%
    filter(max(C08leaf_0h, C08leaf_1h, C08root_0h, C08root_1h)>=5)

pas <- pas %>%
    mutate(C08leaf_0h_cpm = C08leaf_0h/sum(C08leaf_0h)*1000000,
           C08leaf_1h_cpm = C08leaf_1h/sum(C08leaf_1h)*1000000,
           C08root_0h_cpm = C08root_0h/sum(C08root_0h)*1000000,
           C08root_1h_cpm = C08root_1h/sum(C08root_1h)*1000000) %>%
    mutate(mean_cpm = (C08leaf_0h_cpm + C08leaf_1h_cpm + C08root_0h_cpm + C08root_1h_cpm)/4)

## Update function to select two PASs with largest distance
firstTwoPAS <- function(PASlist, count){
    ## PASlist was order by coordinates
    strand = PASlist %>% dplyr::select(strand) %>% unique()
    if(count >=2){
        if(strand=="+"){
            return(bind_rows(PASlist[count,],PASlist[1,]))
        }else{
            return(bind_rows(PASlist[1,],PASlist[count,]))
        }
    }else{
        return(PASlist[1,])
    }
}

pas2 <- pas %>%
    arrange(geneID, start) %>%
    group_nest(geneID) %>%
    mutate(PAS_count = map_int(data, function(x) nrow(x))) %>%
    mutate(data2 = map2(.x=data, .y=PAS_count, .f=firstTwoPAS)) %>%
    select(geneID, PAS_count, data2) %>%
    unnest(data2)



write.table(pas2, file = "selected_PAS_cpm.tsv",
            sep = "\t", quote = F, row.names = F,
            col.names = T)
