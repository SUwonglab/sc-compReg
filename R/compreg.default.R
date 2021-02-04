compreg.default <- function(k,
                            peak.gene.prior.path,
                            thresh,
                            ...) {
    if (is(symbol, 'list')) {
        symbol <- unlist(symbol, use.names=F)
    }

    symbol = unlist(symbol, use.names = F)
    a <- maxk(tf.binding, 5000, 2)

    tf.binding[tf.binding - a[, ncol(a)] < 0] <- 0

    file <- compRegLoad(peak.gene.prior.path)
    f <- match(file$C1, elem.name, nomatch=0)
    d <- f > 0
    f1 <- match(file$C2, symbol, nomatch=0)
    d1 <- f1 > 0
    indices <- (d * d1) == 1
    fc <- cbind(f[indices], f1[indices])

    f2 <- unique(fc, orient = 'r')
    ic <- match(do.call(paste, data.frame(fc)), do.call(paste, data.frame(f2)), nomatch=0)
    ic <- as.numeric(ic)

    c3 <- accumArrayMin(ic, as.numeric(file$C3[indices]))
    c4 <- accumArrayMin(ic, as.numeric(file$C4[indices]))

    c4[c4 < thresh] <- 0
    d0 <- 500000
    c <- exp(-1 * c3 / d0) * c4
    beta <- sparseMatrix(dims = c(length(symbol),length(elem.name)),
                         i = f2[, 2],
                         j = f2[, 1],
                         x = c)

}
