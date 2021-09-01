save.image("logs/build_library.RData")

log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")


library(BSgenome.Hsapiens.NCBI.GRCh38)
library(Biostrings)
library(plyranges)
library(magrittr)
library(tidyverse)

source("lib/helpers.R")

set.seed(snakemake@config$seed)

n_sub_libraries <- snakemake@config$sub_libraries

max_oligo_length <- snakemake@config$oligo$max_len

re_seqs <- DNAStringSet(unlist(snakemake@config$restriction_enzymes))

primer_5p <- snakemake@config$oligo$primer_5p
primer_3p <- snakemake@config$oligo$primer_3p

cloning_site_1 <- re_seqs[snakemake@config$cloning_site$enzyme1]
cloning_site_2 <- re_seqs[snakemake@config$cloning_site$enzyme2]

stuffer_bp <- snakemake@config$cloning_site$stuffer_bp
stuffer_gc <- snakemake@config$cloning_site$stuffer_GC

cloning_site_length <- str_length(cloning_site_1) + stuffer_bp + str_length(cloning_site_2)

barcodes_per_frag <- snakemake@config$bc_per_frag

bc_length <- snakemake@config$oligo$bc_len

max_indel_size <- snakemake@config$oligo$max_indel

random_control_n <- snakemake@config$random_controls
random_control_w_mut_n <- snakemake@config$random_controls_with_mutation

wiggle_room <- max_indel_size + 1

final_frag_len <- max_oligo_length - str_length(primer_5p) - cloning_site_length - str_length(primer_3p) - bc_length - wiggle_room

window_flank <- floor((final_frag_len - 1) / 2)

snps <- read_tsv(snakemake@input$filtered_snps)

snps_cleaned <- snps %>%
    filter(!is.na(ref), !is.na(alt),
           abs(str_count(ref, "[ACGT]") - str_count(alt, "[ACGT]")) <= max_indel_size) %>%
    # mutate(alt = str_extract(alt, "^[ACGT-]+"), # only use first alternate allele if multiple listed
    mutate(alt = str_split(alt, ",")) %>%
    unnest(alt) %>%
    mutate(is_indel = str_count(ref, "[ACGT]") != str_count(alt, "[ACGT]"),  # annotate indels
           is_indel_fixed_base = is_indel & !str_detect(ref, "-") & !str_detect(alt, "-"), # annotate indels with the preceding fixed base
           ref = ifelse(is_indel_fixed_base, str_sub(ref, 2), str_replace(ref, "-", "")), # remove fixed first base or the "-" annotation for indels
           alt = ifelse(is_indel_fixed_base, str_sub(alt, 2), str_replace(alt, "-", ""))) %>%
    mutate(pos = ifelse(is_indel, as.numeric(pos) + 1, as.numeric(pos))) %>%
    distinct(snp, chr, pos, ref, alt, .keep_all = T) %>%
    group_by(snp, chr, pos, ref) %>%
    mutate(alt_n = as.integer(table(alt)[alt])) %>%
    slice_max(alt_n, n = 1, with_ties = F) %>%
    ungroup() %>%
    select(-alt_n) %>%
    arrange(chr, pos, snp) %>%
    mutate(fragment = row_number()) # account for removing fixed base


snps_gr <- snps_cleaned %>%
    transmute(seqnames = chr, start = pos, end = pos + str_length(ref) - 1,
              fragment, chr, pos, snp, ref, alt, is_indel, is_indel_fixed_base, source) %>%
    as_granges() %>%
    mutate()


snps_windows_gr <- resize(snps_gr, width = 2 * window_flank + 1, fix = "center")

## Get a window of sequence around SNP

snp_fragments_raw <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, `seqlevelsStyle<-`(snps_windows_gr, "NCBI"))


snps_frags <- snps_cleaned %>%
    select(fragment, chr, pos, snp, ref, alt, is_indel, is_indel_fixed_base, source) %>%
    mutate(seq = as.character(snp_fragments_raw),
           window_start = start(snps_windows_gr),
           window_end = end(snps_windows_gr))



# remove 10bp from each flank
ignore_edge_bps <- 10
overlapping_snps_gr <- resize(snps_windows_gr, width = width(snps_windows_gr) - 2*ignore_edge_bps,
                              fix = "center") %>%
    select(fragment) %>%
    join_overlap_left(select(snps_gr, -fragment))



# interim_data <- "./outs/LDhap"
# dir.create(interim_data)
# haplo_results_file <- "./outs/haplo_results.rds"


# if (!file.exists(haplo_results_file)) {

#     haplotypes_df <- overlapping_snps_gr %>%
#         as_tibble() %>%
#         filter(!is.na(snp), str_detect(snp, "rs\\d+")) %>%
#         # filter(!is.na(rsID_ldlink), str_detect(rsID_ldlink, "rs\\d+")) %>%
#         group_by(fragment) %>%
#         filter(dplyr::n() > 1) %>%
#         nest() %>%
#         mutate(out_file = file.path(interim_data, paste0("fragment_", fragment, ".txt")),
#                query = map(data, ~ pull(., snp)))
#                # query = map(data, ~ ifelse(str_detect(.x$rsID_ldlink, "rs\\d+"), .x$rsID_ldlink, paste(.x$chr, .x$pos, sep = ":"))))

#     haplo_results <- haplotypes_df %>%
#         mutate(LDhap_data = map2(query, out_file, query_ldhap))


#     write_rds(haplo_results, haplo_results_file)
# } else {
#     haplo_results <- read_rds(haplo_results_file)
# }


# tidy_LDhap_result <- function(x) {
#     haplotype_allele_df <- x %>% mutate(haplotype = row_number()) %>%
#         dplyr::rename(count = Count, frequency = Frequency) %>%
#         gather("snp", "allele", -haplotype, -count, -frequency) %>%
#         group_by(snp, allele) %>%
#         mutate(allele_count = sum(count)) %>%
#         group_by(snp) %>%
#         mutate(allele_freq = allele_count / sum(allele_count[!duplicated(allele)]))
#     return(haplotype_allele_df)
# }


# haplo_results_tidy <- haplo_results %>%
#     mutate(LDhap_data_tidy = map(LDhap_data, tidy_LDhap_result)) %>%
#     select(fragment, LDhap_data_tidy) %>%
#     unnest(LDhap_data_tidy) %>%
#     group_by(fragment) %>%
#     filter(n_distinct(snp) > 1) %>%
#     ungroup()


# haplo_fragments <-
#     haplo_results_tidy %>%
#     mutate(allele = str_replace(allele, "-", "")) %>%
#     left_join(select(snps_cleaned, snp, chr, pos, ref, alt, is_indel)) %>%
#     left_join(select(snps_frags, fragment, seq, window_start, window_end)) %>%
#     filter(allele == ref | allele == alt) %>%
#     mutate(seq_sub = str_sub(seq, pos - window_start + 1, pos - window_start + str_length(ref))) %>%
#     group_by(fragment, haplotype, frequency, count, seq) %>% nest() %>%
#     mutate(haplotype_seq = map2_chr(seq, data, ~ mutate_fragments(.x, .y$pos-.y$window_start+1, .y$ref, .y$allele))) %>%
#     ungroup()

# haplo_fragments_check <- haplo_fragments %>% unnest(data) %>% ungroup() %>%
#     group_by(fragment, haplotype) %>%
#     summarise(full_seq_diff = str_length(haplotype_seq[1]) - str_length(seq[1]),
#               var_seq_diff = sum(str_length(allele) - str_length(ref)))

# # This should only return haplotype fragments where SNP positions are overlapping
# haplo_fragments_check %>% filter(full_seq_diff != var_seq_diff)






snps_frags_allele <- snps_frags %>% mutate(allele = map2(ref, alt, ~ c(.x, .y))) %>% unnest(allele) %>%
    mutate(var_seq = pmap_chr(.,
                              function(seq, pos, window_start, ref, allele, ...) {
                                  mutate_fragments(seq, pos - window_start + 1, ref, allele)
                              })) %>%
    mutate(variant = ifelse(allele == ref, "ref", "alt"))



flank_width <- max(width(re_seqs)) - 1 # for checking RE sites


combined_fragments <- bind_rows(
    snps_frags_allele %>%
        mutate(haplotype = NA) %>%
        select(fragment, variant, haplotype, seq, var_seq, snp, chr, pos, allele, ref, alt, is_indel, window_start, window_end) %>%
        group_by(fragment, variant, haplotype, seq, var_seq) %>% nest() %>% ungroup #,
    # haplo_fragments %>%
    #     filter(haplotype_seq %!in% snps_frags_allele$var_seq,
    #            frequency >= 0.01) %>%
    #     unnest(data) %>%
    #     select(fragment, haplotype, seq, var_seq = haplotype_seq, snp, chr, pos, allele, ref, alt, is_indel, window_start, window_end) %>%
    #     group_by(fragment, haplotype, seq, var_seq) %>% nest() %>% ungroup
) %>%
    mutate(var_seq_flanks  = paste0(str_sub(primer_5p, -flank_width),
                                    var_seq,
                                    str_sub(cloning_site_1, 1, flank_width)))



# Randomly mutate restriction sites in ref/alt pairs
# If haplotype created unique restriction site, toss it


## Search from R.E. sites
re_matches <- combined_fragments %>%
    ungroup %>%
    mutate(re_matches = map(var_seq_flanks, function(x) {
        map_dfr(names(re_seqs), function(re_name) {
            vmatchPattern(re_seqs[[re_name]], x)[[1]] %>%
                as.data.frame() %>% mutate(RE = re_name)
        })
    }))



# re_matches %>% ungroup %>%
#     select(fragment, variant, haplotype, var_seq, re_matches) %>%
#     unnest(re_matches, keep_empty = T) %>%
#     group_by(fragment) %>%
#     filter(any(!is.na(RE))) %>%
#     filter(any(is.na(RE))) %>% View




combined_fragments_corrected <- re_matches %>% # filter(Fragment == 340) %>%
    group_by(fragment) %>%
    group_modify(~ mutate_restriction_sites(.x, .y$fragment))

tmp <- combined_fragments_corrected %>% filter(!fail) %>%
    mutate(var_seq_corrected_flanks  = paste0(str_sub(primer_5p, - flank_width),
                                              var_seq_corrected,
                                              str_sub(cloning_site_1, 1, flank_width))) %>%
    ungroup %>%
    mutate(re_matches = map(var_seq_corrected_flanks, function(x) {
        map_dfr(names(re_seqs), function(re_name) {
            vmatchPattern(re_seqs[[re_name]], x)[[1]] %>%
                as.data.frame() %>% mutate(RE = re_name)
        })
    }))


fragments_final <- combined_fragments_corrected %>% ungroup() %>%
    filter(!fail) %>%
    select(fragment_initial = fragment, variant, haplotype_initial = haplotype, seq, var_seq, var_seq_corrected, data) %>%
    mutate(fragment = dense_rank(fragment_initial)) %>%
    group_by(fragment) %>%
    mutate(haplotype = dense_rank(haplotype_initial)) %>%
    ungroup() %>%
    mutate(frag_id = ifelse(!is.na(variant),
                            paste0('fragment-', str_pad(fragment, width = 4, pad = 0), "_", "allele-", variant),
                            paste0('fragment-', str_pad(fragment, width = 4, pad = 0), "_", "haplotype-", haplotype)))


random_controls <- fragments_final %>% filter(!is.na(variant)) %>%
    distinct(seq) %>%
    # sample_n(random_control_n*2, replace = random_control_n*2 > nrow(.)) %>%
    mutate(fragment_initial = row_number(),
           GC = str_count(seq, "[GC]")/str_length(seq),
           shuffle_seq = map_chr(seq, ~ str_split(., "") %>% unlist %>% sample() %>% str_c(collapse = "")),
           shuffle_GC = str_count(shuffle_seq, "[GC]")/str_length(shuffle_seq),
           pos = floor((str_length(shuffle_seq) + 1)/2),
           ref = str_sub(shuffle_seq, pos, pos),
           alt = ref) %>%
    mutate(variant = map2("ref", "alt", ~ c(.x, .y))) %>% unnest(variant) %>%
    mutate(allele = ifelse(variant == "ref", ref, alt)) %>%
    mutate(var_seq = shuffle_seq) %>%
    mutate(var_seq_flanks  = paste0(str_sub(primer_5p, - flank_width),
                                    var_seq,
                                    str_sub(cloning_site_1, 1, flank_width))) %>%
    ungroup %>%
    mutate(re_matches = map_lgl(var_seq_flanks, function(x) {
        map_dfr(names(re_seqs), function(re_name) {
            vmatchPattern(re_seqs[[re_name]], x)[[1]] %>%
                as.data.frame() %>% mutate(RE = re_name)
        }) %>% nrow %>% is_greater_than(0)
    })) %>% group_by(seq) %>%
    filter(!any(re_matches)) %>%
    ungroup()


random_controls_w_mut <- fragments_final %>% filter(!is.na(variant)) %>%
    distinct(seq) %>%
    # sample_n(random_control_w_mut_n*2, replace = random_control_w_mut_n*2 > nrow(.)) %>%
    mutate(fragment_initial = row_number(),
           GC = str_count(seq, "[GC]")/str_length(seq),
           shuffle_seq = map_chr(seq, ~ str_split(., "") %>% unlist %>% sample() %>% str_c(collapse = "")),
           shuffle_GC = str_count(shuffle_seq, "[GC]")/str_length(shuffle_seq),
           pos = floor((str_length(shuffle_seq) + 1)/2),
           ref = str_sub(shuffle_seq, pos, pos),
           alt = map2_chr(shuffle_seq, ref, ~ str_split(.x, "") %>% unlist %>% str_subset(.y, negate = T) %>% sample(1))) %>%
    mutate(variant = map2("ref", "alt", ~ c(.x, .y))) %>% unnest(variant) %>%
    mutate(allele = ifelse(variant == "ref", ref, alt)) %>%
    mutate(var_seq = pmap_chr(.,
                              function(shuffle_seq, pos, ref, allele, ...) {
                                  mutate_fragments(shuffle_seq, pos, ref, allele)
                              })) %>%
    mutate(var_seq_flanks  = paste0(str_sub(primer_5p, - flank_width),
                                    var_seq,
                                    str_sub(cloning_site_1, 1, flank_width))) %>%
    ungroup %>%
    mutate(re_matches = map_lgl(var_seq_flanks, function(x) {
        map_dfr(names(re_seqs), function(re_name) {
            vmatchPattern(re_seqs[[re_name]], x)[[1]] %>%
                as.data.frame() %>% mutate(RE = re_name)
        }) %>% nrow %>% is_greater_than(0)
    })) %>% group_by(seq) %>%
    filter(!any(re_matches)) %>%
    ungroup()


random_controls_sample <- random_controls %>%
    filter(fragment_initial %in% sample(unique(fragment_initial), random_control_n)) %>%
    mutate(fragment = dense_rank(fragment_initial)) %>%
    mutate(frag_id = paste0('randomfragment-', str_pad(fragment, width = 4, pad = 0), "_", "allele-", variant))

random_controls_w_mut_sample <- random_controls_w_mut %>%
    filter(fragment_initial %in% sample(unique(fragment_initial), random_control_w_mut_n)) %>%
    mutate(fragment = dense_rank(fragment_initial)) %>%
    mutate(frag_id = paste0('randomfragmentmut-', str_pad(fragment, width = 4, pad = 0), "_", "allele-", variant))

random_controls_allele  <- bind_rows(random_controls_sample, random_controls_w_mut_sample)

total_frags <- nrow(fragments_final) + nrow(random_controls_allele)

barcode_df <- read_tsv(snakemake@config$barcodes, col_names = F) %>%
    set_colnames("barcode") %>%
    mutate(barcode_w_flanks = paste0(str_sub(cloning_site_2, start = -flank_width),
                                     barcode,
                                     str_sub(primer_3p, start = 1, end = flank_width))) %>%
    bind_cols(map_dfc(as.list(re_seqs), vcountPattern, subject = .$barcode_w_flanks)) %>%
    mutate(match = purrr::reduce(select(., names(re_seqs)), `+`)) %>%
    filter(match == 0)


barcode_lib <- barcode_df %>%
    sample_n(total_frags * barcodes_per_frag) %>% pull(barcode) %>%
    split(rep(1:total_frags, each = 10)) %>%
    map(sort)




barcoded_library <-
    bind_rows(fragments_final %>% select(frag_id, frag_seq = var_seq_corrected),
              random_controls_allele %>% select(frag_id, frag_seq = var_seq)) %>%
    mutate(barcode = barcode_lib) %>%
    unnest(barcode) %>%
    group_by(frag_id) %>%
    mutate(barcode_id = row_number()) %>%
    ungroup

if (stuffer_bp > 0) {
    gc_bp <- round(stuffer_gc * stuffer_bp)
    at_bp <- stuffer_bp - gc_bp
    stuffer_df <- tibble(
        stuffer = map_chr(1:nrow(barcoded_library),
        ~ paste0(sample(c(sample(c("G", "C"), gc_bp, replace = T),
                            sample(c("A", "T"), at_bp, replace = T)), replace = F), collapse = ""))) %>%
        mutate(stuffer_w_flanks = paste0(str_sub(cloning_site_1, start = 2),
                                     stuffer,
                                     str_sub(cloning_site_2, end = -2)))  %>%
        bind_cols(map_dfc(as.list(re_seqs), vcountPattern, subject = .$stuffer_w_flanks)) %>%
        mutate(match = purrr::reduce(select(., names(re_seqs)), `+`))


    while(any(stuffer_df$match > 0)) {
        stuffer_df_match <- stuffer_df %>% filter(match > 0)

        stuffer_df_match <- tibble(
            stuffer = map_chr(1:nrow(stuffer_df_match),
                              ~ paste0(sample(c(sample(c("G", "C"), gc_bp, replace = T),
                                                sample(c("A", "T"), at_bp, replace = T)), replace = F), collapse = ""))) %>%
            mutate(stuffer_w_flanks = paste0(str_sub(cloning_site_1, start = 2),
                                             stuffer,
                                             str_sub(cloning_site_2, end = -2)))  %>%
            bind_cols(map_dfc(as.list(re_seqs), vcountPattern, subject = .$stuffer_w_flanks)) %>%
            mutate(match = purrr::reduce(select(., names(re_seqs)), `+`))

        stuffer_df <- bind_rows(stuffer_df %>% filter(match == 0),
                                stuffer_df_match)
    }

    stuffer_seqs <- stuffer_df$stuffer

} else {
    stuffer_seqs <- ""
}


final_library <- barcoded_library %>%
    mutate(stuffer = stuffer_seqs) %>%
    mutate(oligo_id = paste0(frag_id, "_barcode-", str_pad(barcode_id, width = 2, pad = 0)),
           oligo = paste0(primer_5p, frag_seq, cloning_site_1, stuffer, cloning_site_2, barcode, primer_3p))



if (n_sub_libraries > 1) {

    final_library <- final_library %>%
        mutate(frag_pair_id = str_extract(frag_id, "^fragment-\\d+"))

    sub_libraries <- final_library %>%
        filter(!is.na(frag_pair_id)) %>%
        distinct(frag_pair_id) %>%
        mutate(sub_library = sample(row_number() %% n_sub_libraries) + 1) %>%
        left_join(final_library %>% filter(!is.na(frag_pair_id)))

    random_controls_sub <- final_library %>%
        filter(str_detect(frag_id, "^randomfragment"))

    sub_libraries <- sub_libraries %>%
        group_by(sub_library) %>% group_nest() %>%
        mutate(data = map(data, ~ bind_rows(., random_controls_sub)))


    walk2(sub_libraries$data, sub_libraries$sub_library,
          ~ .x %>% select(oligo_id, oligo) %>%
              write_csv(str_replace(snakemake@output$oligos, "\\.csv",
                                    paste0("_sublib-", .y, ".csv"))))


    walk2(sub_libraries$data, sub_libraries$sub_library,
          ~ .x %>% select(oligo_id, barcode) %>%
              write_csv(str_replace(snakemake@output$barcode_ref, "\\.csv",
                                    paste0("_sublib-", .y, ".csv"))))

}

final_library %>% select(oligo_id, oligo) %>%
    write_csv(snakemake@output$oligo)

final_library %>% select(oligo_id, barcode) %>%
    write_csv(snakemake@output$barcode_ref)






variant_final <- fragments_final %>% filter(!is.na(variant)) %>% unnest(data) %>%
    select(fragment, variant, snp, chr, pos, allele, ref, alt,
           ref_fragment = seq, var_fragment = var_seq_corrected,
           frag_start = window_start, frag_end = window_end, frag_id)
write_csv(variant_final, snakemake@output$variant_ref)

# haplotype_final <- fragments_final %>% filter(!is.na(haplotype)) %>% unnest(data) %>%
#     select(fragment, haplotype, snp, chr, pos, allele, ref, alt,
#            ref_fragment = seq, var_fragment = var_seq_corrected,
#            frag_start = window_start, frag_end = window_end, frag_id)
# write_csv(haplotype_final, "./data/raw/lib3_design/skin_disease_haplotype_ref.csv")

random_controls_final <-  random_controls_allele %>%
    select(fragment, variant, pos, allele, ref, alt, ref_fragment = shuffle_seq, var_fragment = var_seq, frag_id)
write_csv(random_controls_final, snakemake@output$random_ctrl_ref)






