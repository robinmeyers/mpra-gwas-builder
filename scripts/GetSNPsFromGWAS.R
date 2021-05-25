
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

disease_list <-  read_tsv(snakemake@config$disease_list)

gwas_catalog <- read_tsv(snakemake@config$gwas_catalog)

if (snakemake@config$extra_gwas != "") {
    extra_gwas <- read_tsv(snakemake@config$extra_gwas)    
} else {
    extra_gwas <- tibble()
}


gwas_index_snps <- gwas_catalog %>%
    filter(`DISEASE/TRAIT` %in% disease_list$GWAS_term,
           as.numeric(`P-VALUE`) <= as.numeric(snakemake@config$gwas_pvalue_threshold)) %>%
    filter(str_detect(SNPS, "rs\\d+ x rs\\d+", negate = T)) %>%
    bind_rows(extra_gwas) %>%
    left_join(disease_list, by = c("DISEASE/TRAIT" = "GWAS_term"))


gwas_snps_rsID <- gwas_index_snps %>%
    filter(str_detect(SNPS, "^rs\\d+$"))

gwas_snps_no_rsID <- gwas_index_snps %>%
    filter(!str_detect(SNPS, "^rs\\d+$"))


if (nrow(gwas_snps_no_rsID) > 0) {


    gwas_snps_no_rsID_fix <- gwas_snps_no_rsID %>%
        mutate(SNPS_fix = SNPS) %>%
        mutate(SNPS_fix = ifelse(str_detect(SNPS_fix, "^chr[0-9XY]+[^0-9XY]+\\d+$"),
                                 str_match(SNPS_fix, "^(chr[0-9XY]+)[^0-9XY]+(\\d+)$") %>% {paste0(.[,2], ":", .[,3])}, SNPS_fix)) %>%
        mutate(SNPS_fix = ifelse(str_detect(SNPS_fix, "^[0-9XY]+-\\d+$"),
                                 str_match(SNPS_fix, "^([0-9XY]+)-(\\d+)$") %>% {paste0("chr", .[,2], ":", .[,3])}, SNPS_fix))

    ## Manually inspect - anything not in "chr:pos" format will be filtered
    # gwas_snps_no_rsID_fix %>% select(SNPS, SNPS_fix) %>% View

    gwas_snps_no_rsID_fix_done <- gwas_snps_no_rsID_fix %>%
        filter(str_detect(SNPS_fix, "^chr[0-9XY]+:\\d+$")) %>%
        mutate(SNPS = SNPS_fix) %>%
        select(-SNPS_fix)

    gwas_snps_final <- bind_rows(gwas_snps_rsID,
                                 gwas_snps_no_rsID_fix_done)
} else {
    gwas_snps_final <- gwas_snps_rsID
}


## Make sure all SNPS are in either rsID or chr:pos format
all(str_detect(gwas_snps_final$SNPS, "^rs\\d+$") |
        str_detect(gwas_snps_final$SNPS, "^chr[0-9XY]+:\\d+$"))


gwas_snps_final %>% write_tsv(snakemake@output$gwas)


# gwas_snps_final %>%
#     distinct(SNPS) %>%
#     write_tsv(snakemake@output$index_snps)


