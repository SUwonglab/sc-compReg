subpopulation_link.default <- function(E.mean.healthy,
                                       E.mean.cll,
                                       O.mean.healthy,
                                       O.mean.cll) {
    if (! is(E.mean.healthy, 'matrix') & ! is(E.mean.healthy, 'Matrix')) {
        stop('E.mean.healthy must be a matrix. Please check your format.')
    }

    if (! is(E.mean.cll, 'matrix') & ! is(E.mean.cll, 'Matrix')) {
        stop('E.mean.cll must be a matrix. Please check your format.')
    }

    if (! is(O.mean.healthy, 'matrix') & ! is(O.mean.healthy, 'Matrix')) {
        stop('O.mean.healthy must be a matrix. Please check your format.')
    }

    if (! is(O.mean.cll, 'matrix') & ! is(O.mean.cll, 'Matrix')) {
        stop('O.mean.cll must be a matrix. Please check your format.')
    }

    this.call <- match.call()

    r1 <- cor(E.mean.healthy, E.mean.cll)
    r2 <- cor(O.mean.healthy, O.mean.cll)
    rr1 <- r1 - t(colSums(r1) %*% t(rowSums(r1))) / sum(r1)
    rr2 <- r2 - t(colSums(r2) %*% t(rowSums(r2))) / sum(r2)
    rr <- rr1 + rr2

    b <- which(rr > 0)
    a <- matrix(0, length(b),7)
    a[, 1:2] <- which(rr > 0, arr.ind = T)
    a[, 3] <- r1[b]
    a[, 4] <- r2[b]
    a[, 5] <- rr1[b]
    a[, 6] <- rr2[b]
    a[, 7] <- rr[b]
    f <- order(a[, 7], decreasing = T)
    a <- a[f, ]
    a <- a[(((a[, 5] > 0) * (a[, 6] > 0)) == 1),  ]

    a.dim <- dim(a)[1]
    match.mat <- matrix(NA, a.dim, 7)
    S1 <- matrix(NA, a.dim, 1)
    S2 <- matrix(NA, a.dim, 1)
    match.idx <- 1
    for (i in 1:a.dim) {
        idx <- is.element(a[i, 1], S1) + is.element(a[i, 2], S2)
        S1[i] <- a[i, 1]
        S2[i] <- a[i, 2]
        if (idx < 2) {
            match.mat[match.idx, ] <- a[i, ]
            match.idx = match.idx + 1
        }
    }
    match.mat <- match.mat[rowSums(is.na(match.mat)) != ncol(match.mat), ]

    output <- list()
    output$match <- match.mat
    output$call <- this.call

    return(output)
}
