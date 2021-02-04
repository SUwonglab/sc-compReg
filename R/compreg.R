compreg <- function(k,
                    peak.gene.prior.path,
                    thresh,
                    p.val.thresh = 0.1) {
    UseMethod("compreg.default")
}
