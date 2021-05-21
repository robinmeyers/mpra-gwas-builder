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



if ("epigenome_csv" %in% names(snakemake@config) && file.exists(snakemake@config$epigenome_csv)) {

    epigenome_csv <- read_csv(snakemake@config$epigenome_csv)
    epigenome_filter_keys <- epigenome_csv %>%
    	filter(filter) %>% pull(name)

} else {

    epigenome_keys <- names(snakemake@config$epigenome)
    epigenome_filter_keys <- epigenome_keys[map(snakemake@config$epigenome, ~ .$filter)]
}

print(epigenome_filter_keys)

ld_snps_filter <- snps_epigenome %>%
    filter(map_lgl(str_split(Epigenome, "\\|"), ~ any(. %in% epigenome_filter_keys))) %>%  # |
               # (Epigenome == "H3K27me3" & !is.na(eQTL))) %>%
    arrange(chr, pos)

write_tsv(ld_snps_filter, snakemake@output$filtered_snps)

