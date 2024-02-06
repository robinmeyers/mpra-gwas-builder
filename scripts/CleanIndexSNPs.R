## Long-running script due to many API queries

if (!interactive()) {
    readr::write_rds(snakemake, paste0(snakemake@log[[1]], ".snakemake.rds"))
    log <- file(snakemake@log[[1]], open="wt")
    sink(log, type = "message")
    sink(log, type = "output")
} 


library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)
library(tidyverse)

source("lib/helpers.R")
source("lib/get_pops.R")
source("lib/harmonize_snps.R")


hg19_to_hg38_chain <- import.chain("assets/hg19ToHg38.over.chain")

# threads <- 4

# if (threads > 1) {
#     library(doMC)
#     registerDoMC(cores = threads)
#     do_parallel <- T
# } else {
#     do_parallel <- F
# }

index_snp_table <- read_tsv(snakemake@input$gwas,
                            col_types = cols(.default = col_character()), quote = "")
# index_snps <- read_tsv("./data/raw/lib3_design/skin_disease_index_snps.txt")



# if (!is.null(snakemake@config$extra_gwas) && snakemake@config$extra_gwas != "") {
#     extra_gwas <- read_tsv(snakemake@config$extra_gwas,
#         col_types = cols(.default = col_character()), quote = "")
# } else {
#     extra_gwas <- tibble()
# }

all(str_detect(index_snp_table$SNPS, "^rs\\d+$") |
        str_detect(index_snp_table$SNPS, "^chr[0-9XY]+:\\d+$"))


index_snps <- index_snp_table %>%
    select(disease = Disease, gwas_snp = SNPS, chr = CHR_ID, pos = CHR_POS,
           pubmed = PUBMEDID, sample = `INITIAL SAMPLE SIZE`) %>%
    mutate(coord = ifelse(is.na(chr), NA, paste0("chr", chr, ":", pos))) %>%
    mutate(coord = ifelse(is.na(coord) & str_detect(gwas_snp, "chr.+:\\d+"), gwas_snp, coord))

index_snps_gr <- index_snps %>%
    filter(!is.na(coord)) %>%
    extract(coord, c("chr", "pos"), "chr([0-9XY]+):([0-9]+)") %>%
    mutate(start = pos, end = pos) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
# index_snps_gr <- index_snps %>%
#     filter(!is.na(coord_b38)) %>%
#     extract(coord_b38, c("chr", "pos"), "chr([0-9XY]+):([0-9]+)") %>%
#     mutate(start = pos, end = pos) %>%
#     makeGRangesFromDataFrame(keep.extra.columns = T)

snps_find_rsid_b37 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, index_snps_gr)
snps_find_rsid_b38 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38, index_snps_gr)
snps_find_rsid_b38_xtra <- snpsByOverlaps(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38,
                                          `seqlevelsStyle<-`(index_snps_gr, "dbSNP")) %>%
    `seqlevelsStyle<-`("NCBI")



snps_find_rsid_b37_tbl <-
   as.data.frame(snps_find_rsid_b37) %>%
    mutate(coord_b37 = paste0("chr", seqnames, ":", pos)) %>%
    select(rs_id_rescue_b37 = RefSNP_id, coord_b37)

snps_find_rsid_b38_tbl <-
    bind_rows(as.data.frame(snps_find_rsid_b38) %>%
                  mutate(coord_b38 = paste0("chr", seqnames, ":", pos)),
              as.data.frame(snps_find_rsid_b38_xtra) %>%
                  mutate(coord_b38 = paste0("chr", seqnames, ":", start))) %>%
    select(rs_id_rescue = RefSNP_id, coord_b38)
# snps_find_rsid_b38_tbl <- snps_find_rsid_b38 %>% as.data.frame() %>%
#     mutate(coord_b38 = paste0("chr", seqnames, ":", pos)) %>%
#     select(rs_id_rescue = RefSNP_id, coord_b38)

index_snps_harmonized <- left_join(index_snps, snps_find_rsid_b38_tbl, by = c("coord" = "coord_b38"), keep = T) %>%
    left_join(snps_find_rsid_b37_tbl, by = c("coord" = "coord_b37"), keep = T) %>%
    mutate(pmap_dfr(list(gwas_snp, coord, rs_id_rescue, coord_b38, rs_id_rescue_b37, coord_b37),
                   harmonize_snps)) %>%
    select(disease, gwas_snp, index_snp, coord_b38, coord_b37, pubmed, sample)


if (!is.null(snakemake@config$gwas_pop_key)) {
    
    study_key_table <-
        index_snps_harmonized %>%
            distinct(pubmed, sample) %>% 
            get_pops_from_samples(snakemake@config$gwas_pop_key)

    index_snps_pop_match <- index_snps_harmonized %>%
        left_join(study_key_table) %>%
        distinct() %>%
        group_by(disease, gwas_snp, index_snp, coord_b38, coord_b37, pubmed, sample) %>%
        summarise(pops = paste0(sort(unique(unlist(str_split(code, ",")))), collapse = ",")) %>%
        ungroup()


    # write_tsv(index_snps_pop_match, "outs/gwas_study_index_snps_matched_populations.tsv")

    index_snps_pop_match %>%
        group_by(disease, pubmed, sample, pops) %>%
        summarise(n_snps = n_distinct(index_snp, na.rm = T)) %>%
        write_tsv("outs/gwas_study_matched_populations.tsv")


} else {
    index_snps_pop_match <- index_snps_harmonized %>%
        mutate(pops = NA_character_)
    
    # index_snps_pop_match <- tibble(disease = character(),
    #                                pubmed = character(),
    #                                sample = character(),
    #                                index_snp = character(),
    #                                pops = character())
}


if (!is.null(snakemake@config$extra_gwas) && snakemake@config$extra_gwas != "") {
    extra_gwas <- read_tsv(snakemake@config$extra_gwas,
        col_types = cols(.default = col_character())) %>%
        mutate(pops = str_replace_all(pops, " ", ""))
    

} else {
    extra_gwas <- tibble()
}

index_snps_extra <-
    bind_rows(index_snps_pop_match, extra_gwas)

rs_ids <- index_snps_extra %>%
    filter(is.na(coord_b38)) %>%
    pull(index_snp) %>%
    unique()

extra_snps_b38 <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38,
    rs_ids, ifnotfound = "drop")
extra_snps_b38xtra <- snpsById(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38,
    rs_ids, ifnotfound = "drop") %>%
    `seqlevelsStyle<-`("NCBI")

index_snps_cleaned <- index_snps_extra %>%
    left_join(as.data.frame(extra_snps_b38) %>%
                    transmute(index_snp = RefSNP_id,
                            coord_b38_rescue = paste0("chr", seqnames, ":", pos)),
                by = "index_snp") %>%
    left_join(as.data.frame(extra_snps_b38xtra) %>%
                    transmute(index_snp = RefSNP_id,
                            coord_b38extra_rescue = paste0("chr", seqnames, ":", start)),
                by = "index_snp") %>%
    mutate(coord_b38 = coalesce(coord_b38, coord_b38_rescue, coord_b38extra_rescue))

index_snps_unique <- index_snps_cleaned %>%
    group_by(index_snp) %>%
    summarise(pops = sort(unique(unlist(str_split(pops, ",")))) %>%
        str_subset(pattern = "^$", negate = T) %>%
        paste0(collapse = ","))

index_snps_unique %>%
    write_csv(snakemake@output$index_snps_unique)

index_snps_cleaned %>%
    select(disease, gwas_snp, index_snp, coord_b38, coord_b37, pubmed, pops, sample) %>%
    write_csv(snakemake@output$index_snps)

