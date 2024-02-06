get_pops_from_samples <- function(study_samples, gwas_pop_key_file) {
    gwas_pop_key <- read_tsv(gwas_pop_key_file)
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


    sample_terms <- study_samples %>%
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
    return(study_key_table)
}