## Long-running script due to many API queries


save.image("logs/clean_index_snps.RData")

log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")


library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)
library(tidyverse)

source("lib/helpers.R")


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

all(str_detect(index_snp_table$SNPS, "^rs\\d+$") |
        str_detect(index_snp_table$SNPS, "^chr[0-9XY]+:\\d+$"))


index_snps <- index_snp_table %>%
    select(disease = Disease, gwas_snp = SNPS, chr = CHR_ID, pos = CHR_POS,
           pubmed = PUBMEDID, sample = `INITIAL SAMPLE SIZE`) %>%
    mutate(coord_b38 = ifelse(is.na(chr), NA, paste0("chr", chr, ":", pos))) %>%
    mutate(coord_b38 = ifelse(is.na(coord_b38) & str_detect(gwas_snp, "chr.+:\\d+"), gwas_snp, coord_b38))


index_snps_gr <- index_snps %>%
    filter(!is.na(coord_b38)) %>%
    extract(coord_b38, c("chr", "pos"), "chr([0-9XY]+):([0-9]+)") %>%
    mutate(start = pos, end = pos) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

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

index_snps_cleaned <- left_join(index_snps, snps_find_rsid_b38_tbl) %>%
    mutate(index_snp = ifelse(str_detect(gwas_snp, "^rs\\d+"), gwas_snp,
                              ifelse(!is.na(rs_id_rescue), rs_id_rescue, NA))) %>%
    left_join(snps_find_rsid_b37_tbl, by = c("coord_b38" = "coord_b37")) %>%
    mutate(coord_b37 = ifelse(is.na(index_snp) & !is.na(rs_id_rescue_b37), coord_b38, NA),
           coord_b38 = ifelse(is.na(index_snp) & !is.na(rs_id_rescue_b37), NA, coord_b38),
           index_snp = ifelse(is.na(index_snp) & !is.na(rs_id_rescue_b37), rs_id_rescue_b37, index_snp)) %>%
    mutate(index_snp = ifelse(is.na(index_snp), gwas_snp, index_snp)) %>%
    select(disease, gwas_snp, index_snp, coord_b38, coord_b37, pubmed, sample)


write_csv(index_snps_cleaned, snakemake@output$index_snps)

