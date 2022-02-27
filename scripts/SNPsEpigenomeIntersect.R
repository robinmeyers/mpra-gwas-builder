
save.image("logs/intersect_epigenome.RData")


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

set.seed(snakemake@config$seed)

## Load data
# ldlink_full_results <- read_tsv("./data/raw/lib3_design/ldlink_full_results.txt")
# haploreg_full_results <- read_tsv("./data/raw/lib3_design/haploreg_full_results.txt")

ld_snps <- read_tsv(snakemake@input$ld_snps)

hg19_to_hg38_chain <- import.chain("assets/hg19ToHg38.over.chain")

if ("epigenome_csv" %in% names(snakemake@config) && file.exists(snakemake@config$epigenome_csv)) {

    epigenome_csv <- read_csv(snakemake@config$epigenome_csv)
    epigenome_keys <- epigenome_csv$name
    epigenome_bed <- map2(epigenome_csv$bedfile, epigenome_csv$genome, function(bedfile, genome) {
        bed <- read_narrowpeaks(bedfile)
        if (genome == "hg19") {
            bed <- liftOver(bed, hg19_to_hg38_chain) %>% unlist
        }
        return(bed)
        })

} else {

    epigenome_keys <- names(snakemake@config$epigenome)
    epigenome_bed <- map(snakemake@config$epigenome, function(epigenome) {
        bed <- read_narrowpeaks(epigenome$bedfile)
        if (epigenome$genome == "hg19") {
            bed <- liftOver(bed, hg19_to_hg38_chain) %>% unlist
        }
        return(bed)
        })
}

epigenome_df <- tibble(key = epigenome_keys,
                       bed = epigenome_bed) %>%
    mutate(key = str_replace_all(key, "[^A-Za-z0-9_]", "_")) %>%
    group_by(key) %>%
    summarise(bed = list(reduce(bed, union_ranges)))

epigenome_keys <- epigenome_df$key
epigenome_bed <- epigenome_df$bed %>% set_names(epigenome_df$key)

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


if (!is.null(snakemake@config$eqtls)) {
    eqtls <-
        map_dfr(snakemake@config$eqtls,
            ~ read_tsv(.$file), .id = "tissue")

    eqtls <- eqtls %>%
        extract(variant_id, c("chr", "pos"), "^(chr[0-9XY]+)_(\\d+)", remove = F) %>%
        mutate(pos = as.integer(pos))
} else {
    eqtls <- tibble(chr = character(), pos = integer(), variant_id = character())
}




ld_snps_epigenome <- ld_snps_gr %>%
    as_tibble(.name_repair = "minimal") %>%
    select(-end, -width, -strand) %>%
    dplyr::rename(chr = seqnames,
                  pos = start) %>%
    mutate(across(all_of(epigenome_keys), ~ ifelse(!is.na(.), cur_column(), NA), .names = "{.col}_dummy")) %>%
    unite(Epigenome, ends_with("_dummy"), sep = ";", na.rm = T) %>%
    left_join(eqtls %>% distinct(chr, pos, eQTL = variant_id))


write_tsv(ld_snps_epigenome, snakemake@output$epigenome)


peak_stats <- tibble(peakset = epigenome_keys, bed = epigenome_bed) %>%
    mutate(peak_num = map_int(epigenome_bed, length),
           peak_width = map_int(epigenome_bed, ~ sum(width(.)))) %>%
    select(-bed)

write_csv(peak_stats, snakemake@output$peak_stats)




