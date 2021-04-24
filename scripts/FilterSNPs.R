log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

# TODO: generalize the filtering process, allow user to specify filters


# index_snps <- read_tsv("./data/raw/lib3_design/skin_disease_gwas.tsv")
#
# snps_ldlink <- read_tsv("./data/raw/lib3_design/SNPS_LDlink/full_results.txt")
# snps_haploreg <- read_tsv("./data/raw/lib3_design/SNPS_HaploReg/full_results.txt")

snps_epigenome <- read_tsv(snakemake@input$epigenome)



ld_snps_filter <- snps_epigenome %>%
    filter(str_detect(Epigenome, "(ATAC|H3K27ac|H3K4me1)")) %>%  # |
               # (Epigenome == "H3K27me3" & !is.na(eQTL))) %>%
    arrange(chr, pos)

write_tsv(ld_snps_filter, snakemake@output$filtered_snps)

