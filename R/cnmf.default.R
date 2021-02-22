cnmf.default <- function(PeakO,
                         X,
                         D,
                         k,
                         beta_min,
                         alpha=0.5,
                         beta_max_scale=5,
                         verbose=T) {
    if (! is(PeakO, 'sparseMatrix')) {
        PeakO <- Matrix(PeakO, sparse = T)
    }

    if (! is(X, 'sparseMatrix')) {
        X <- Matrix(X, sparse = T)
    }

    if (! is(D, 'sparseMatrix')) {
        D <- Matrix(D, sparse = T)
    }

    if (! is(k, 'numeric')) {
        stop('k must be an integer.')
    }

    k <- as.integer(k)

    if (!missing(beta_min)) {
        if (! is(beta_min, numeric)) {
            warning('User must supply a valid beta_min of type numeric, or supply nothing at all. Using default value.')
            beta <- 10.0^seq(beta_max_scale, -3, -1)
        } else {
            beta <- beta_min * 10.0^seq(beta_max_scale, 0, -1)
        }
    } else {
        beta <- 10.0^seq(beta_max_scale, -3, -1)
    }
    # convert to matrix to use NNLM::nnmf function
    X.nmf <- nnmf(as.matrix(X), k=k)
    w2 <- X.nmf$W
    h2 <- X.nmf$H

    PO.nmf <- nnmf(as.matrix(PeakO), k=k)
    w1 <- PO.nmf$W
    h1 <- PO.nmf$H

    PO.dim <- dim(PeakO)

    mat.init <- initializeMatrix(PO.dim[1],
                          PO.dim[2],
                          dim(X)[2],
                          k,
                          D)

    for (j in 1:length(beta)) {
        dp.ret <- compute_lambda(PeakO,
                                w1,
                                h1,
                                X,
                                w2,
                                h2,
                                D,
                                alpha,
                                beta[j])
        lambda1 <- dp.ret$lambda1
        lambda2 <- dp.ret$lambda2
        if (verbose) {
            message(paste("On ", j, "-th loop", sep = ''))
            message(paste("Lambda 1:", lambda1))
            message(paste("Lambda 2:", lambda2))
        }
        nmf.ret <- nmf_cluster_sep2_lap(PeakO,
                             X,
                             D,
                             k,
                             max_iter=300,
                             lambda1,
                             lambda2,
                             mat.init$W1,
                             mat.init$H1,
                             mat.init$W2,
                             mat.init$H2,
                             verbose=verbose)
        update.ret <- postLapMatMult(nmf.ret$W1,
                       nmf.ret$W2,
                       nmf.ret$H1,
                       nmf.ret$H2)
        H1 <- update.ret$H1
        H2 <- update.ret$H2
        W1 <- update.ret$W1
        W2 <- update.ret$W2
        score <- nmf.ret$score
    }
    cluster.output <- cluster(H1, H2);
    return(list("W1" = W1,
                "W2" = W2,
                "H1" = H1,
                "H2" = H2,
                "lambda1" = lambda1,
                "lambda2" = lambda2,
                "score" = score,
                "cluster1" = cluster.output$c1,
                "cluster2" = cluster.output$c2))
}
