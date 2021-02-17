cluster.profile <- function(O1,
                            E1,
                            O1.idx,
                            E1.idx,
                            symb1,
                            pk.name1,
                            O2,
                            E2,
                            O2.idx,
                            E2.idx,
                            symb2,
                            pk.name2,
                            pk.name.intersect1,
                            pk.name.intersect2,
                            ...) {
    UseMethod("cluster.profile.default")
}
