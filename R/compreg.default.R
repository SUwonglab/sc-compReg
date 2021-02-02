compreg.default <- function(k,
                            peak.gene.prior.path) {
    if (is(symbol, 'list')) {
        symbol <- unlist(symbol, use.names=F)
    }

    a <- maxk(tf.binding, k, 2) # applied over each row
    tf.binding[tf.binding - a[, k] < 0] <- 0
    file <- compRegLoad(peak.gene.prior.path)
    f <- match(file$C1, elem.name, nomatch=0)
    d <- f > 0
    f1 <- match(file$C2, symbol)
    d1 <- f1 > 0

}
