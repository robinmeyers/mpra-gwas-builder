save.image("logs/filter_snps.RData")

log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

set.seed(snakemake@config$seed)

snps_epigenome <- read_tsv(snakemake@input$epigenome)

if (!is.null(snakemake@config$txdb_filters)) {
    txdb_filter_keys <- snakemake@config$txdb_filters
} else {
    txdb_filter_keys <- character()
}


if ("epigenome_csv" %in% names(snakemake@config) && file.exists(snakemake@config$epigenome_csv)) {

    epigenome_csv <- read_csv(snakemake@config$epigenome_csv)
    epigenome_filter_keys <- epigenome_csv %>%
    	filter(filter) %>% pull(name) %>%
        str_replace_all("[^A-Za-z0-9_]", "_")

} else {

    epigenome_keys <- names(snakemake@config$epigenome)
    epigenome_filter_keys <- epigenome_keys[map_lgl(snakemake@config$epigenome, ~ .$filter)] %>%
        str_replace_all("[^A-Za-z0-9_]", "_")
}



print(txdb_filter_keys)
print(epigenome_filter_keys)

ld_snps_filter <- snps_epigenome %>%
    filter(map_lgl(str_split(Epigenome, ";"), ~ any(. %in% epigenome_filter_keys))) %>%
    filter(map_lgl(str_split(txdb_annot, ";"), ~ ! any(. %in% txdb_filter_keys))) %>%
    arrange(chr, pos)

write_tsv(ld_snps_filter, snakemake@output$filtered_snps)

