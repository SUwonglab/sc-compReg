compreg <- function(symbol,
                    tf.binding,
                    elem.name,
                    peak.gene.prior.path,
                    E1,
                    E1.idx,
                    E2,
                    E2.idx,
                    O1.mean,
                    O2.mean,
                    match.mat,
                    thresh = 0.2,
                    sig.level = 0.05,
                    top.tf = 5000,
                    d0.default = 500000,
                    ...) {
    UseMethod("compreg.default")
}
