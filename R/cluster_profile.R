"cluster_profile" <- function(O1,
                              E1,
                              O1.idx,
                              E1.idx,
                              symbol1,
                              peak.name1,
                              O2,
                              E2,
                              O2.idx,
                              E2.idx,
                              symbol2,
                              peak.name2,
                              peak.name.intersect1,
                              peak.name.intersect2,
                              ...) {
    UseMethod("cluster_profile")
}
