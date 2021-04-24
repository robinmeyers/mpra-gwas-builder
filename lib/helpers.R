`%!in%` <- Negate(`%in%`)

# liftover_range <- function(range, chain) {
#     tibble(range = range) %>%
#         extract(range, c("chr", "start", "end"), "^(chr[0-9XY]+):(\\d+)-(\\d+)") %>%
#         makeGRangesFromDataFrame() %>%
#         liftOver(chain) %>%
#         unlist %>%
#         as.data.frame() %>%
#         mutate(range = paste0(seqnames, ":", start, "-", end)) %>%
#         pull(range)
# }
#
# liftover_position <- function(coord, chain) {
#     pos <- tibble(coord = coord) %>%
#         extract(coord, c("chr", "pos"), "^(chr[0-9XY]+):(\\d+)", remove = F)
#
#     pos_filtered <- pos %>% filter(!is.na(coord))
#     pos_filtered$coord_new <-
#         GPos(pos_filtered$chr, pos_filtered$pos, coord = pos_filtered$coord) %>%
#         liftOver(chain) %>%
#         unlist %>% as("GPos") %>%
#         as.data.frame() %>%
#         mutate(coord = paste0(seqnames, ":", pos)) %>%
#         pull(coord)
#
#     pos %>% left_join(pos_filtered, by = "coord") %>% pull(coord_new)
# }

query_ldlink <- function(snp, pop, out_dir, r2 = 0.8, force = F, retry_errors = T) {
    out_file <- file.path(out_dir, paste0(pop, "_",
                                          str_replace(snp, ":", "-"), ".txt"))


    if (file.exists(out_file) && !force) {
        results <- read_tsv(out_file, col_types = cols())

        if ("RS_Number" %in% colnames(results)) {
            results %>% 
                mutate(RegulomeDB = as.character(RegulomeDB)) %>%
                return()

        }

        if (!retry_errors) {
            return(results %>% mutate(RegulomeDB = NA_character_))
        }
            
    }

    raw_results <- LDproxy(snp = snp, pop = pop, token = Sys.getenv("LDLINK_TOKEN"))

    if ("R2" %!in% colnames(raw_results)) {
        error_results <- raw_results %>%
            mutate(SNP = snp, Population = pop,
                   Error = T)
        write_tsv(error_results, out_file)
        return(error_results)
    }

    filtered_results <- raw_results %>%
        filter(R2 >= r2) %>%
        mutate(RegulomeDB = as.character(RegulomeDB))

    write_tsv(filtered_results, out_file)
    return(filtered_results)
}

query_haploreg <- function(snps, pop, out_dir, r2 = 0.8, force = F) {
    cat(pop, "\n")
    out_file <- file.path(out_dir, paste0(pop, "_haploreg.tsv"))
    
    if (file.exists(out_file) && !force) {
        haploreg_results <- read_tsv(out_file, col_types = cols())
        return(haploreg_results)
    }


    haploreg_results <-
        queryHaploreg(snps, ldPop = pop, ldThresh = r2,  timeout = 1000000) %>%
        mutate(Population = pop,
               chr = as.character(chr))
    write_tsv(haploreg_results, out_file)
    return(haploreg_results)
}


mutate_fragments <- function(seq, pos, ref, allele) {
    # seq <- "ACGTCGCGTCACTGGCT"
    # pos <- c(4, 8, 15)
    # ref <- c("T", "GT", "")
    # alt <- c("A", "", "A")

    # seq <- test_hap$seq[1] %>% str_sub(85,105)
    # pos <- test_hap$pos - test_hap$window_start + 1 - 84
    # ref <- test_hap$ref
    # allele <- test_hap$Allele

    get_ref <- str_sub(seq, pos, pos + str_length(ref) - 1)

    alt_i <- allele != get_ref

    if (!any(alt_i)) return(seq)

    pos <- pos[alt_i]
    ref <- ref[alt_i]
    allele <- allele[alt_i]

    order_i <- order(pos, str_length(ref))
    pos <- pos[order_i]
    ref <- ref[order_i]
    # alt <- alt[order_i]
    allele <- allele[order_i]

    # print(pos)
    # print(ref)
    # print(allele)

    var_sites <- IRanges(start = pos, end = pos + str_length(ref) - 1)
    var_sites_1 <- var_sites
    start(var_sites_1) <- start(var_sites) - 1
    end(var_sites_1) <- end(var_sites) + 1

    var_bps <- IRanges::reduce(var_sites)
    var_site_which <- findOverlaps(var_bps, var_sites_1)
    var_site_alleles <- split(allele[subjectHits(var_site_which)], queryHits(var_site_which)) %>%
        map_chr(paste0, collapse = "")
    # print(var_bps)

    fixed_starts <- c(1, end(var_bps) + 1)
    fixed_ends <- c(start(var_bps) - 1, str_length(seq))
    fixed_subseqs <- str_sub(seq, fixed_starts, fixed_ends)


    # print(fixed_starts)
    # print(fixed_ends)

    # print(var_site_alleles)

    zipped_alt <- zipup(fixed_subseqs, c(var_site_alleles, "")) %>% unlist %>% paste0(collapse = "")
    #print(pairwiseAlignment(seq, zipped_alt))
    return(zipped_alt)
}


mutate_restriction_sites <- function(x, frag) {
    print(frag)
    # print(re_sites)
    re_sites <- x %>% unnest(re_matches)

    # Condition with no restriction sites
    if (nrow(re_sites) == 0) {
        corrected_x <- x %>%
            mutate(fail = F,
                   var_seq_corrected = var_seq)
        return(corrected_x)
    }

    # print(all(!is.na(re_sites$Haplotype)))
    # Condition where the only RE sites are in haplotype fragments
    if (all(!is.na(re_sites$haplotype))) {
        corrected_x <- x %>% mutate(fail = !is.na(haplotype) & (haplotype %in% re_sites$haplotype),
                                    var_seq_corrected = ifelse(fail, NA, var_seq))
        return(corrected_x)
    }

    # Condition where ref or alt creates/breaks RE site
    if (re_sites %>% filter(!is.na(variant)) %>%
        count(variant, RE) %>% complete(variant = c("ref", "alt"), RE, fill = list(n = 0)) %>%
        spread(variant, n) %$% any(alt != ref)) {
        corrected_x <- x %>% mutate(fail = T, var_seq_corrected = NA)
        return(corrected_x)
    }

    # Condition where there's at least one of ref and alt with a restriction site
    var_df <- x %>% unnest(data)


    alignments <- map2(x$seq, x$var_seq, pairwiseAlignment)
    aligned_seq <-  map(alignments, alignedPattern)
    aligned_var_seq <- map(alignments, alignedSubject)
    align_map_df <- map(alignments,
                        ~tibble(query = alignedSubject(.) %>% as.matrix() %>% str_detect("[ACGT]") %>% cumsum %>% add(flank_width),
                                ref = alignedPattern(.) %>% as.matrix() %>% str_detect("[ACGT]") %>% cumsum))

    align_q_r_maps <- map(align_map_df, ~ distinct(., query, .keep_all = T) %>% select(query, ref) %>% deframe)
    align_r_q_maps <- map(align_map_df, ~ distinct(., ref, .keep_all = T) %>% select(ref, query) %>% deframe)


    re_sites <- x %>%
        mutate(re_matches = map2(re_matches, align_q_r_maps,
                                 ~ mutate(.x,
                                          r_start = .y[as.character(start)],
                                          r_end = .y[as.character(end)]) %>%
                                     mutate(r_start = ifelse(is.na(r_start), 1, r_start),
                                            r_end = ifelse(is.na(r_end), str_length(x$seq[1]), r_end)))) %>%
        unnest(re_matches)

    re_sites_iranges <- IRanges(re_sites$r_start, re_sites$r_end)
    var_iranges <- IRanges(var_df$pos - var_df$window_start + 1, var_df$pos - var_df$window_start + str_length(var_df$ref))

    seq <- x$seq[1]
    var_seq <- x$seq

    re_sites_dist <- distanceToNearest(re_sites_iranges, var_iranges)
    mcols(re_sites_iranges)$distToVar <- mcols(re_sites_dist)$distance

    proposed_subs <- tibble(pos = numeric(), ref = character(), alt = character())

    # start with RE site closest to any variant pos
    for (re_site_i in order(mcols(re_sites_iranges)$distToVar)) {

        new_seq <- seq %>% str_split("") %>% unlist() %>%
            magrittr::inset(proposed_subs$pos, proposed_subs$alt) %>% str_c(collapse = "")

        new_subseq <- str_sub(new_seq, re_sites$r_start[re_site_i], re_sites$r_end[re_site_i])
        if (re_sites$r_start[re_site_i] == 1 & str_length(new_subseq) < str_length(re_seqs[re_sites$RE[re_site_i]])) {
            new_subseq <- str_c(str_sub(primer_5p, -(str_length(re_seqs[re_sites$RE[re_site_i]]) - str_length(new_subseq))),
                                new_subseq)
        }
        if (re_sites$r_end[re_site_i] == str_length(new_seq) & str_length(new_subseq) < str_length(re_seqs[re_sites$RE[re_site_i]])) {
            new_subseq <- str_c(new_subseq,
                                str_sub(cloning_site, 1, str_length(re_seqs[re_sites$RE[re_site_i]]) - str_length(new_subseq)))
        }

        # check if the site still exists with the proposed subs
        if (new_subseq != re_seqs[re_sites$RE[re_site_i]]) {
            next
        }

        # Nominate furthest bp to mutate
        re_site <- re_sites_iranges[re_site_i]
        re_pos_iranges <- IRanges(start(re_site):end(re_site), start(re_site):end(re_site))
        re_pos_iranges_dist <- distanceToNearest(re_pos_iranges, var_iranges)
        re_pos_i <- start(re_pos_iranges)[which.max(mcols(re_pos_iranges_dist)$distance)]
        re_pos_bp <- str_sub(new_seq, re_pos_i, re_pos_i)


        for (new_bp in sample(setdiff(c("A", "C", "G", "T"), re_pos_bp))) {
            proposed_new_seq <- str_c(str_sub(new_seq, re_pos_i - flank_width, re_pos_i - 1),
                                      new_bp,
                                      str_sub(new_seq, re_pos_i + 1, re_pos_i + flank_width))
            if (re_pos_i - flank_width < 1) {
                proposed_new_seq <- str_c(str_sub(primer_5p, -(flank_width - re_pos_i + 1)), proposed_new_seq)
            }
            if (re_pos_i + flank_width > str_length(new_seq)) {
                proposed_new_seq <- str_c(proposed_new_seq, str_sub(cloning_site, 1, flank_width + re_pos_i - str_length(new_seq)))
            }

            if (any(map_lgl(as.character(re_seqs), ~ str_detect(proposed_new_seq, .)))) {
                # try next bp if the first didn't work
                next
            } else {
                # if it did work list the new subs and move to next RE site
                proposed_subs <- proposed_subs %>%
                    bind_rows(tibble(pos = re_pos_i, ref = re_pos_bp, alt = new_bp))
                break
            }
            # else fail the fragment
            corrected_x <- x %>% mutate(fail = T, var_seq_corrected = fail)
            return(corrected_x)

        }

    }


    ## Make the proposed subs to the reference sequence and each variant
    corrected_ref_seq <- seq %>% str_split("") %>% unlist() %>%
        magrittr::inset(proposed_subs$pos, proposed_subs$alt) %>% str_c(collapse = "")
    ## Double check none of them contain re sites



    final_test_RE_sites <- x %>%
        mutate(var_seq_corrected = map2_chr(corrected_ref_seq, data,
                                            ~ mutate_fragments(.x, .y$pos-.y$window_start+1, .y$ref, .y$allele))) %>%
        mutate(contains_RE = map_lgl(var_seq_corrected, ~ any(str_detect(., as.character(re_seqs)))))

    ## If ref or alt do, fail the entire fragment, otherwise just fail the haplotype
    if (nrow(final_test_RE_sites %>% filter(!is.na(variant), contains_RE)) > 0) {
        corrected_x <- final_test_RE_sites %>%
            mutate(var_seq_corrected = NA,
                   fail = T) %>%
            select(-contains_RE)
    } else {
        corrected_x <- final_test_RE_sites %>%
            mutate(var_seq_corrected = ifelse(contains_RE, NA, var_seq_corrected),
                   fail = contains_RE) %>%
            select(-contains_RE)
    }

    return(corrected_x)

}
