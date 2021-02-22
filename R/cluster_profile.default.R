cluster_profile.default <- function(O1,
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
                                    peak.name.intersect2) {
    this.call <- match.call()

    K1 <- max(max(O1.idx), max(E1.idx))
    symbol <- intersect(symbol1, symbol2)
    f1 <- match(symbol, symbol1, nomatch=0)
    E1.mean <- matrix(0, length(f1), K1)
    for (i in 1:K1) {
        gp <- which(E1.idx == i)
        if (length(gp) == 0) next
        E1.mean[, i] <- Matrix::rowMeans(E1[f1, gp])
    }

    K2 <- max(max(O2.idx), max(E2.idx))
    f2 <- match(symbol, symbol2, nomatch=0)
    E2.mean <- matrix(0, length(f2), K2)
    for (i in 1:K2) {
        gp = which(E2.idx == i)
        if (length(gp) == 0) next
        E2.mean[, i] <- Matrix::rowMeans(E2[f2, gp])
    }

    pf1 <- match(peak.name.intersect1, peak.name1, nomatch=0)
    pf2 <- match(peak.name.intersect2, peak.name2, nomatch=0)
    d1 <- is.element(peak.name1, peak.name.intersect1)
    d2 <- is.element(peak.name2, peak.name.intersect2)
    elem.name <- c(peak.name.intersect1, peak.name1[d1 == 0],
                  peak.name2[d2 == 0])

    m <- length(peak.name.intersect1)
    elem.len <- length(elem.name)
    O1.mean <- matrix(0, elem.len, K1)
    O2.mean <- matrix(0, elem.len, K2)

    sum.d1 <- sum(d1 == 0)
    for (i in 1:K1) {
        gp = which(O1.idx == i)
        if (length(gp) == 0) next
        O1.mean[1:m, i] = Matrix::rowMeans(O1[pf1, gp] > 0)
        O1.mean[(1+m) : (m + sum.d1), i] = Matrix::rowMeans(O1[d1 == 0, gp] > 0)
    }

    for (i in 1:K2) {
        gp = which(O2.idx == i)
        if (length(gp) == 0) next
        O2.mean[1:m, i] = Matrix::rowMeans(O2[pf2, gp] > 0)
        O2.mean[(m + sum.d1 + 1) : elem.len, i] = Matrix::rowMeans(O2[d2 == 0, gp] > 0)
    }

    O1.mean = sweep(O1.mean, 2, Matrix::colMeans(O1.mean), '/')
    O2.mean = sweep(O2.mean, 2, Matrix::colMeans(O2.mean), '/')

    output <- list('E1.mean' = E1.mean,
                   'E2.mean' = E2.mean,
                   'symbol' = symbol,
                   'O1.mean' = O1.mean,
                   'O2.mean' = O2.mean,
                   'elem.name' = elem.name)
    output$call <- this.call

    return(output)
}
