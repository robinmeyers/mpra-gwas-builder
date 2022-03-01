
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)
library(cowplot)

disease_list <- read_csv(snakemake@config$disease_list)

gwas_snps_raw <- read_tsv(snakemake@input$gwas_snps)
index_snps_raw <- read_csv(snakemake@input$index_snps, col_types = 'ccccc')
ld_snps_raw <- read_tsv(snakemake@input$ld_snps)
epigenome_snps_raw <- read_tsv(snakemake@input$epigenome_snps)
filtered_snps_raw <- read_tsv(snakemake@input$filtered_snps)
final_snps_raw <- read_csv(snakemake@input$variant_ref)

fig_dir <- "outs/figures"
dir.create(fig_dir, recursive = T)


gwas_snps <- full_join(gwas_snps_raw, index_snps_raw, by = c("Disease" = "disease", "SNPS" = "gwas_snp"))

gwas_study_count <- gwas_snps %>% distinct(Disease, PUBMEDID) %>% count(Disease)
gwas_study_count %>%
    ggplot(aes(x = fct_rev(Disease), y = n)) +
    geom_col() +
    coord_flip() +
    labs(y = "GWAS Studies",
         title = paste(n_distinct(gwas_snps$PUBMEDID, na.rm = T), "GWAS Studies")) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
ggsave(file.path(fig_dir, "gwas_studies.pdf"))

gwas_snp_count <- gwas_snps %>% distinct(Disease, index_snp) %>%
    group_by(Disease) %>%
    summarise(n = n_distinct(index_snp, na.rm = T))
gwas_snp_count %>%
    ggplot(aes(x = fct_rev(Disease), y = n)) +
    geom_col() +
    coord_flip() +
    labs(y = "GWAS SNPs",
         title = paste(n_distinct(gwas_snps$SNPS, na.rm = T), "GWAS Index SNPs")) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
ggsave(file.path(fig_dir, "index_snps_per_disease.pdf"))

linked_snps <- gwas_snps %>%
    select(Disease, index_snp) %>%
    full_join(ld_snps_raw, by = "index_snp")

linked_snp_count <- linked_snps %>% distinct(Disease, snp) %>%
    group_by(Disease) %>%
    summarise(n = n_distinct(snp, na.rm = T))
linked_snp_count %>%
    ggplot(aes(x = fct_rev(Disease), y = n)) +
    geom_col() +
    coord_flip() +
    labs(y = "Linked SNPs",
         title = paste(n_distinct(linked_snps$snp, na.rm = T), "SNPs in LD (R2 > 0.8) with Index SNPs")) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
ggsave(file.path(fig_dir, "linked_snps_per_disease.pdf"))



linked_snps %>%
    distinct(index_snp, snp) %>%
    group_by(index_snp) %>%
    summarise(n = n_distinct(snp, na.rm = T)) %>%
    ggplot(aes(x = n)) +
    geom_histogram() +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    labs(x = "Number of linked SNPs",
         y = "Index SNPs",
         title = "Distribution of linked SNPs per index SNP") +
    theme_cowplot()
ggsave(file.path(fig_dir, "linked_snps_per_index_snp.pdf"))



filtered_snps <- gwas_snps %>%
    select(Disease, index_snp) %>%
    full_join(filtered_snps_raw, by = "index_snp")

filtered_snp_count <- filtered_snps %>% distinct(Disease, snp) %>%
    group_by(Disease) %>%
    summarise(n = n_distinct(snp, na.rm = T))
filtered_snp_count %>%
    ggplot(aes(x = fct_rev(Disease), y = n)) +
    geom_col() +
    coord_flip() +
    labs(y = "Filtered SNPs",
         title = paste(n_distinct(filtered_snps$snp, na.rm = T), "SNPs post filtering")) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
ggsave(file.path(fig_dir, "filtered_snps_per_disease.pdf"))


filtered_snps %>% distinct(index_snp, snp) %>%
    group_by(index_snp) %>%
    summarise(n = n_distinct(snp, na.rm = T)) %>%
    ggplot(aes(x = n)) +
    geom_histogram() +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    labs(x = "Number of filtered SNPs",
         y = "Index SNPs",
         title = "Distribution of filtered SNPs per index SNP") +
    theme_cowplot()
ggsave(file.path(fig_dir, "filtered_snps_per_index_snp.pdf"))


epigenome_snps <- epigenome_snps_raw %>%
    select(snp, Epigenome) %>%
    mutate(Epigenome = str_split(Epigenome, ";")) %>%
    unnest(Epigenome) %>%
    filter(!is.na(Epigenome)) %>%
    distinct()

epigenome_snps %>%
    ggplot(aes(x = fct_rev(Epigenome))) +
    geom_bar() +
    coord_flip() +
    labs(x = "Peak set", y = "Filtered SNPs") +
    theme_cowplot()
ggsave(file.path(fig_dir, "filtered_snps_per_peakset.pdf"))


peak_stats <- read_csv(snakemake@input$peak_stats)
epigenome_stats <- epigenome_snps %>%
    count(Epigenome) %>%
    left_join(peak_stats, by = c("Epigenome" = "peakset")) %>%
    mutate(snps_per_peak = n / peak_num,
           snps_per_mb = n / peak_width * 1e6)

epigenome_stats %>%
    ggplot(aes(x = fct_rev(Epigenome), y = snps_per_mb)) +
    geom_col() +
    coord_flip() +
    labs(x = "Peak set", y = "SNPs per Mb") +
    theme_cowplot()
ggsave(file.path(fig_dir, "filtered_snps_per_peakset_mb.pdf"))

final_snps <- filtered_snps %>%
    select(Disease, index_snp, snp) %>%
    full_join(final_snps_raw)

final_snp_count <- final_snps %>% distinct(Disease, fragment) %>%
    group_by(Disease) %>%
    summarise(n = n_distinct(fragment, na.rm = T))


final_snp_count %>%
    ggplot(aes(x = fct_rev(Disease), y = n)) +
    geom_col() +
    coord_flip() +
    labs(y = "Final SNPs",
         title = paste(n_distinct(final_snps$fragment, na.rm = T), "SNPs in final library")) +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
ggsave(file.path(fig_dir, "final_snps_per_disease.pdf"))


final_snps %>%
    distinct(index_snp, fragment) %>%
    group_by(index_snp) %>%
    summarise(n = n_distinct(fragment, na.rm = T)) %>%
    ggplot(aes(x = n)) +
    geom_histogram() +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    labs(x = "Number of linked SNPs",
         y = "Index SNPs",
         title = "Distribution of linked SNPs per index SNP in final library") +
    theme_cowplot()
ggsave(file.path(fig_dir, "final_snps_per_index_snp.pdf"))


stats_table <- gwas_snp_count %>% rename(GWAS_SNPs = n) %>%
    left_join(linked_snp_count %>% rename(Linked_SNPs = n)) %>%
    left_join(filtered_snp_count %>% rename(Filtered_SNPs = n)) %>%
    left_join(final_snp_count %>% rename(Final_SNPs = n))


stats_table_unique <- tibble(Disease = "Unique") %>%
    mutate(GWAS_SNPs = n_distinct(gwas_snps$index_snp, na.rm=T),
           Linked_SNPs = n_distinct(linked_snps$snp, na.rm=T),
           Filtered_SNPs = n_distinct(filtered_snps$snp, na.rm=T),
           Final_SNPs = n_distinct(final_snps$fragment, na.rm=T))

bind_rows(stats_table, stats_table_unique) %>%
    write_csv(snakemake@output$stats)

