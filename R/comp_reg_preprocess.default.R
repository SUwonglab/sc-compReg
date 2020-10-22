comp_reg_preprocess.default <- function(peak.name.intersect.dir,
                                        token='\t') {
    # if (class(sample1) == "character") {
    #     # load sample1 from directory
    #     sample1 = readMat(sample1)
    # }
    #
    # if (class(sample2) == "character") {
    #     # load sample2 from directory
    #     sample2 = readMat(sample2)
    # }
    file = load_pn_intersect_file(peak.name.intersect.dir, token)
    peak.name.intersect = cbind(file$vo, file$vt)
    return(peak.name.intersect)
}
