compreg <- function(k,
                    peak.gene.prior.path,
                    thresh = 0.2,
                    p.val.thresh = 0.1,
                    sig.level = 0.05,
                    top.tf = 5000,
                    d0.default = 500000) {
    UseMethod("compreg.default")
}
