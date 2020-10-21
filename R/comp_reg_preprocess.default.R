comp_reg_preprocess.default <- function(E1,
                                        E1.idx,
                                        O1,
                                        O1.idx,
                                        peak.name1,
                                        symbol1,
                                        E2,
                                        E2.idx,
                                        O2,
                                        O2.idx,
                                        peak.name2,
                                        symbol2,
                                        peak.name.intersect.dir,
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
