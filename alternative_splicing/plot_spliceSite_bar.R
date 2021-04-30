#! /usr/bin/env Rscript

library(tidyverse)


d <- read_tsv("splice_sites_stat.input.tsv",
              col_names = c("chr", "start", "end",
                            "group", "Type"))
d <- d %>%
    mutate(Type = factor(Type,
                         levels = c("Annotated",
                                    "ShortReads",
                                    "Novel")))
count_tbl <- d %>%
    group_by(group, Type) %>%
    summarize(count = n())

group_sum <- count_tbl %>%
    group_by(group) %>%
    summarize(total=sum(count))

canonical_t <- group_sum %>%
    filter(group == "Canonical") %>% pull(total)

non_canonical_t <- group_sum %>%
    filter(group == "Non-Canonical") %>% pull(total)

count_tbl <- count_tbl %>%
    mutate(percentage = case_when(
               group == "Canonical" ~ count/canonical_t,
               group == "Non-Canonical" ~ count/non_canonical_t))

print(count_tbl)
write_tsv(count_tbl,
          file = "spliceSite_stat.out")

p <- ggplot(d, aes(x=group,
                   fill=Type))

p <- p + geom_bar(position = "fill",
                  color = "black",
                  width = 0.6,
                  alpha = 0.8)

p <- p + scale_fill_brewer(palette = "Dark2")
p <- p + xlab("") + ylab("Percentage")
p <- p + theme_bw()
p <- p + theme(panel.grid = element_blank(),
               ##legend.position = "top",
               axis.text = element_text(color="black"),
               axis.text.x = element_text(angle= 90,
                                          vjust = 0.5,
                                          hjust =1))

ggsave(filename = "spliceSite_annotated.bar.pdf",
       width = 3.5, height = 3.5)
