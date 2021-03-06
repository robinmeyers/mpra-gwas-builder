## Long-running script due to many API queries

readRenviron(".Renviron")

save.image("logs/get_snps_in_ld.RData")

log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

if (! "haploR" %in% rownames(installed.packages())) {
    options(repos = list(CRAN="http://cran.rstudio.com/"))
    install.packages("haploR")
}

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(LDlinkR)
library(haploR)
library(VariantAnnotation)
library(magrittr)
library(tidyverse)

source("lib/helpers.R")

set.seed(snakemake@config$seed)

hg19_to_hg38_chain <- import.chain("assets/hg19ToHg38.over.chain")

# threads <- 4

# if (threads > 1) {
#     library(doMC)
#     registerDoMC(cores = threads)
#     do_parallel <- T
# } else {
#     do_parallel <- F
# }

index_snps_cleaned <- read_csv(snakemake@input$index_snps)
# index_snps <- read_tsv("./data/raw/lib3_design/skin_disease_index_snps.txt")

r2_threshold <- snakemake@config$r2_threshold
r2_threshold_pop_specific <- snakemake@config$r2_threshold_pop_spec


pops <- snakemake@config$pops
# pops <- c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL")

if (!is.null(snakemake@config$gwas_pop_key)) {
    gwas_pop_key <- read_tsv(snakemake@config$gwas_pop_key)

    sample_types <- c("individuals?",
                      "cases?",
                      "controls?",
                      "men",
                      "women",
                      "boys?",
                      "girls?",
                      "adults?",
                      "adolescents?",
                      "children and adolescents",
                      "children",
                      "infants?",
                      "neonates?",
                      "mothers?",
                      "fathers?",
                      "parents?",
                      "males?",
                      "females?",
                      "users?",
                      "non-users?",
                      "families",
                      "trios?",
                      "responders?",
                      "non-responders?",
                      "attempters?",
                      "nonattempters?",
                      "alcohol drinkers?",
                      "drinkers?",
                      "non-drinkers?",
                      "smokers?",
                      "non-smokers?",
                      "donors?",
                      "twin pairs?",
                      "twins?",
                      "child sibling pairs?",
                      "fetuses",
                      "offspring",
                      "early adolescents?",
                      "remitters?",
                      "non-remitters?",
                      "athletes?",
                      "Individuals?",
                      "indivduals?",
                      "triads?",
                      "patients?",
                      "pairs?",
                      "case-parent trios?",
                      "recipients?",
                      "affected child",
                      "long sleepers?",
                      "short sleepers?",
                      "unaffected relatives?",
                      "carriers?",
                      "non-carriers?",
                      "cell lines?",
                      "indiviudals?",
                      "referents?",
                      "individuuals?",
                      "duos?",
                      "indivdiuals?",
                      "inidividuals?")

    number_regex <- "(?:(?<=(?:\\s|\\b))\\d+(?:\\,\\d+)*(?=\\s))"
    type_regex <- paste0("(?:", paste0(sample_types, collapse = "|"), ")")


    full_regex <- paste0(
        "(", number_regex, ")", # greedy match first number
        "\\s*((?:(?!.*", type_regex, ").*)|(?:.*?))\\s*", # Greedy match rest if no sample type in lookahead, or passive match
        "(", type_regex,  "?(?!.*", type_regex, "))") # Match last sample type by ensuring no sample type in lookahead

    # split_regex <- "(?<!\\d)(,[\\s\\,]*| and )(?=[\\sA-Aa-z]*[0-9]+[,0-9]*[0-9]+\\s)"
    split_regex <- paste0("((?:,+[,\\s]*\\s+)|(?:and ))(?=[\\sA-Aa-z]*", number_regex, ")")


    sample_terms <- index_snps_cleaned %>%
        distinct(pubmed, sample) %>%
        mutate(split_sample = str_split(sample, split_regex)) %>%
        unnest(split_sample)


    full_matches <- bind_cols(sample_terms,
                              str_match(sample_terms$split_sample, full_regex) %>%
                                  set_colnames(c("match", "number", "capture", "type")) %>%
                                  as_tibble())

    study_key_table <- full_matches %>%
        distinct(pubmed, sample, split_sample, capture) %>%
        rename(term = capture) %>%
        left_join(gwas_pop_key) %>%
        filter(!is.na(code))

    index_snps_pop_match <- index_snps_cleaned %>%
        left_join(study_key_table) %>%
        distinct() %>%
        group_by(disease, gwas_snp, index_snp, coord_b38, coord_b37, pubmed, sample) %>%
        summarise(pops = paste0(sort(unique(unlist(str_split(code, ",")))), collapse = ",")) %>%
        ungroup()


    write_tsv(index_snps_pop_match, "outs/gwas_study_index_snps_matched_populations.tsv")

    index_snps_pop_match %>%
        group_by(disease, pubmed, sample, pops) %>%
        summarise(n_snps = n_distinct(index_snp, na.rm = T)) %>%
        write_tsv("outs/gwas_study_matched_populations.tsv")


} else {
    index_snps_pop_match <- tibble(disease = character(),
                                   pubmed = character(),
                                   sample = character(),
                                   index_snp = character(),
                                   pops = character())
}

max_pops <- snakemake@config$max_pops


index_snps_pop_match_filtered <- index_snps_pop_match %>%
    filter(!is.na(pops) & pops != "") %>%
    filter(map_lgl(str_split(pops, ","), ~ length(.) <= max_pops))


index_snps_pop_all <- crossing(index_snp = unique(index_snps_cleaned$index_snp),
                               pop = pops) %>%
    bind_rows(index_snps_pop_match_filtered %>%
                  mutate(pop = str_split(pops, ",")) %>%
                  unnest(pop) %>%
                  distinct(index_snp, pop))


snps_to_query <- index_snps_pop_all %>%
    filter(str_detect(index_snp, "rs\\d+"),
           !is.na(pop) & pop != "") %>%
    mutate(r2_threshold = ifelse(is.null(r2_threshold_pop_specific) | pop == "ALL",
                                 r2_threshold, r2_threshold_pop_specific))

out_dir <- "outs/SNPS_LDlink"

dir.create(out_dir, showWarnings = F, recursive = T)

ldlink_results <- snps_to_query %>%
    mutate(ldlink_results = pmap(list(index_snp, pop, r2_threshold),
        ~ query_ldlink(snp = ..1, pop = ..2, r2 = ..3, out_dir = out_dir, retry_errors = snakemake@config$retry_errors)))

ldlink_results_table <- ldlink_results %>%
    unnest(ldlink_results) %>%
    filter(R2 >= r2_threshold)

write_tsv(ldlink_results_table, "outs/ldlink_full_results.txt")

haploreg_pops <- c("AFR" = "AFR",
                   "AMR" = "AMR",
                   "EAS" = "ASN",
                   "EUR" = "EUR",
                   "SAS" = "ASN")

out_dir_haploreg <- "outs/SNPS_HaploReg"
dir.create(out_dir_haploreg, showWarnings = F, recursive = T)

haploreg_results <- snps_to_query %>%
    filter(pop %in% names(haploreg_pops)) %>%
    mutate(pop = haploreg_pops[pop]) %>%
    group_by(pop, r2_threshold) %>%
    summarise(index_snps = list(sort(index_snp))) %>%
    mutate(haploreg_results = pmap(list(index_snps, pop, r2_threshold),
        ~ query_haploreg(snps = ..1, pop = ..2, r2 = ..3,
                        force = T, out_dir = out_dir_haploreg))) %>%
    ungroup()

if (nrow(haploreg_results) > 0) {
    haploreg_results_table <- haploreg_results %>%
        select(pop, r2_threshold, haploreg_results) %>%
        unnest(haploreg_results) %>%
        select(index_snp = query_snp_rsid, everything()) %>%
        filter(r2 >= r2_threshold)
} else {
    haploreg_results_table <- tibble(
        index_snp = character(),
        pop = character(),
        chr = character(),
        pos_hg38 = character(),
        r2 = double(),
        D = double(),
        is_query_snp = double(),
        rsID = character(),
        ref = character(),
        alt = character()
    )
}

write_tsv(haploreg_results_table, "outs/haploreg_full_results.txt")



# Harmonize rsIDs and genomic coordinates for all LD SNPs from both sources

# ldlink_results_table <- read_tsv("./data/raw/lib3_design/ldlink_full_results.txt")
# haploreg_results_table <- read_tsv("./data/raw/lib3_design/haploreg_full_results.txt")

## LDlink data is in hg19 coordinates
ldlink_snps <- ldlink_results_table %>%
    extract(Alleles, c("ref", "alt"), "([ACGT-]+)\\/([ACGT-]+)", remove = F) %>%
    filter(!is.na(ref), !is.na(alt))

ldlink_snps_b38 <- ldlink_snps %>%
    extract(Coord, c("chr", "start"), "(chr[0-9XY]+):(\\d+)", remove = F) %>%
    mutate(end = start) %>%
    select(seqnames = chr, start, end, snp = RS_Number, index_snp, coord_b37 = Coord, ref, alt) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    liftOver(hg19_to_hg38_chain) %>% unlist %>%
    as_tibble() %>%
    mutate(coord_b38 = paste0(seqnames, ":", start),
           snp = ifelse(is.na(snp) | !str_detect(snp, "^rs\\d+"), coord_b38, snp)) %>%
    select(snp, coord_b38, ref, alt, index_snp, coord_b37) %>%
    distinct()



## HaploReg data is in hg38 coordinates, but not all snps returned have genome coordinates

haploreg_snps <- haploreg_results_table %>%
    mutate(coord_b38 = ifelse(is.na(chr), NA, paste0("chr", chr, ":", pos_hg38))) %>%
    select(snp = rsID, coord_b38, ref, alt, index_snp)

haploreg_snps_no_coord <- haploreg_snps %>%
    filter(is.na(coord_b38)) %>% pull(snp) %>% unique()

## Try to rescue location data from SNPlocs packages and GTEx variant info

haploreg_snps_find_locs_b38 <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, haploreg_snps_no_coord, ifnotfound = "drop") %>%
    GRanges() %>% as_tibble() %>%
    mutate(chr = str_replace(as.character(seqnames), "^(chr|ch)", "")) %>%    select(chr, pos_b38 = start, snp = RefSNP_id)
haploreg_snps_find_locs_b38_xtra <- snpsById(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38, haploreg_snps_no_coord, ifnotfound = "drop") %>%
    GRanges() %>% as_tibble() %>%
    mutate(chr = str_replace(as.character(seqnames), "^(chr|ch)", "")) %>%
    select(chr, pos_b38 = start, snp = RefSNP_id)


if (!is.null(snakemake@config$gtex_table)) {
    gtex_var_map <- read_tsv(snakemake@config$gtex_table,
                             col_types = "c-----cc") %>%
        dplyr::rename(rs_id = "rs_id_dbSNP151_GRCh38p7")

    haploreg_snps_find_locs_gtex <- gtex_var_map %>% filter(rs_id %in% haploreg_snps_no_coord) %>%
        extract(variant_id, c("chr", "pos_b38"), "^chr([0-9XY]+)_(\\d+)") %>%
        mutate(pos_b38 = as.numeric(pos_b38)) %>%
        select(chr, pos_b38, snp = rs_id)
} else {
    haploreg_snps_find_locs_gtex <- tibble()
}


haploreg_snps_find_locs_combined <- bind_rows(
    haploreg_snps_find_locs_b38,
    haploreg_snps_find_locs_b38_xtra,
    haploreg_snps_find_locs_gtex
) %>% distinct %>%
    mutate(coord_b38_rescue = paste0("chr", chr, ":", pos_b38)) %>%
    select(snp, coord_b38_rescue)


haploreg_snps_b38 <- haploreg_snps %>%
    left_join(haploreg_snps_find_locs_combined) %>%
    mutate(coord_b38 = as.character(ifelse(is.na(coord_b38), coord_b38_rescue, coord_b38))) %>%
    select(-coord_b38_rescue) %>%
    distinct()



## Combine LD SNPs


ld_snps_b38 <- bind_rows(
    ldlink_snps_b38 %>% mutate(source = "LDlink"),
    haploreg_snps_b38 %>% mutate(source = "HaploReg")
)



## Get TxDb annotations

ld_snps_b38_gr <- ld_snps_b38 %>%
    extract(coord_b38, c("seqnames", "start"), "(.+):(\\d+)") %>%
    filter(!is.na(seqnames), !is.na(start)) %>%
    mutate(end = start) %>%
    select(seqnames, start, end, snp) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

ld_snps_txdb_loc <- locateVariants(ld_snps_b38_gr, TxDb.Hsapiens.UCSC.hg38.knownGene, AllVariants())

ld_snps_txdb_loc_df <- as_tibble(ld_snps_txdb_loc) %>%
    transmute(coord_b38 = paste0(seqnames, ":", start),
              txdb_annot = LOCATION) %>%
    distinct() %>%
    group_by(coord_b38) %>%
    summarise(txdb_annot = paste0(txdb_annot, collapse = ";"))

ld_snps_b38_annot <- left_join(ld_snps_b38,
                               ld_snps_txdb_loc_df, by = "coord_b38")


write_tsv(ld_snps_b38_annot, snakemake@output$ld_snps)
