log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

## Load packages
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(VariantAnnotation)
# library(rtracklayer)
# library(plyranges)


## Set up project
# library(ProjectTemplate)
# load.project()

# str(snakemake@config$epigenome)

library(rtracklayer)
library(plyranges)
library(tidyverse)

## Load data
# ldlink_full_results <- read_tsv("./data/raw/lib3_design/ldlink_full_results.txt")
# haploreg_full_results <- read_tsv("./data/raw/lib3_design/haploreg_full_results.txt")

ld_snps <- read_tsv(snakemake@input$ld_snps)

hg19_to_hg38_chain <- import.chain("assets/hg19ToHg38.over.chain")

if ("epigenome_csv" %in% names(snakemake@config) && file.exists(snakemake@config$epigenome_csv)) {

    epigenome_csv <- read_csv(snakemake@config$epigenome_csv)
    epigenome_keys <- epigenome_csv$name
    epigenome_bed <- map2(epigenome_csv$bedfile, epigenome_csv$genome, function(bedfile, genome) {
        bed <- read_bed(bedfile)
        if (genome == "hg19") {
            bed <- liftOver(bed, hg19_to_hg38_chain) %>% unlist
        }
        return(bed)
        }) %>% set_names(epigenome_keys)

} else {

    epigenome_keys <- names(snakemake@config$epigenome)
    epigenome_bed <- map(snakemake@config$epigenome, function(epigenome) {
        bed <- read_bed(epigenome$bedfile)
        if (epigenome$genome == "hg19") {
            bed <- liftOver(bed, hg19_to_hg38_chain) %>% unlist
        }
        return(bed)
        })
}




ld_snps_gr <- ld_snps %>%
    filter(!is.na(coord_b38)) %>%
    extract(coord_b38, c("chr", "pos"), "(chr[0-9XY]+):(\\d+)", remove = F) %>%
    mutate(start = pos, end = pos) %>%
    select(-pos) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

epigenome_ranges <- map(epigenome_bed,
    ~ as_tibble(.) %>%
    mutate(range = paste0(seqnames, ":", start, "-", end)) %>%
    pull(range))

mcols(ld_snps_gr) <- cbind(mcols(ld_snps_gr),
    map2_dfc(epigenome_ranges, epigenome_bed, ~ .x[findOverlaps(ld_snps_gr, .y, maxgap = 0, select = "first")]))


eqtls <- 
    map_dfr(snakemake@config$eqtls,
        ~ read_tsv(.$file), .id = "tissue")

eqtls <- eqtls %>%
    extract(variant_id, c("chr", "pos"), "^(chr[0-9XY]+)_(\\d+)", remove = F) %>%
    mutate(pos = as.integer(pos))


ld_snps_epigenome <- ld_snps_gr %>%
    as_tibble() %>%
    select(-end, -width, -strand) %>%
    dplyr::rename(chr = seqnames,
                  pos = start) %>%
    mutate(across(all_of(epigenome_keys), ~ ifelse(!is.na(.), cur_column(), NA), .names = "{.col}_dummy")) %>%
    unite(Epigenome, ends_with("_dummy"), sep = ";", na.rm = T) %>%
    left_join(eqtls %>% distinct(chr, pos, eQTL = variant_id))


write_tsv(ld_snps_epigenome, snakemake@output$epigenome)


