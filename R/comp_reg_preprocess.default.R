comp_reg_preprocess.default <- function(peak.name.intersect.dir,
                                        token='\t') {
    file <- load_pn_intersect_file(peak.name.intersect.dir, token)
    return(file)
}
