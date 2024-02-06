harmonize_snps <- function(snp, coord, rs_id_rescue, coord_b38, rs_id_rescue_37, coord_b37) {
    
    if (!is.na(rs_id_rescue) && snp == rs_id_rescue && coord == coord_b38) {
        return(tibble(index_snp = snp, coord_b38 = coord_b38, coord_b37 = NA))
    } else if (!is.na(rs_id_rescue_37) && snp == rs_id_rescue_37 && coord == coord_b37) {
        return(tibble(index_snp = snp, coord_b38 = NA, coord_b37 = coord_b37))
    } else if (!is.na(coord_b37) && str_detect(snp, "^chr") && coord == coord_b37) {
        return(tibble(index_snp = rs_id_rescue_37, coord_b38 = NA, coord_b37 = coord_b37))
    } else if (!is.na(coord_b38) && coord == coord_b38) {
        return(tibble(index_snp = rs_id_rescue, coord_b38 = coord_b38, coord_b37 = NA))
    } else if (!is.na(coord_b37) && coord == coord_b37) {
        return(tibble(index_snp = rs_id_rescue_37, coord_b38 = NA, coord_b37 = coord_b37))
    } else {
        return(tibble(index_snp = snp, coord_b38 = NA, coord_b37 = NA))
    }
}
