cnmf.default <- function(peak.o,
                         X,
                         D,
                         k,
                         beta.min,
                         rna.cell.label,
                         atac.cell.label,
                         alpha=0.5,
                         beta.max.scale=5,
                         verbose=T,
                         ...) {
    if (! is(peak.o, 'matrix') & ! is(peak.o, 'sparseMatrix')) {
        stop('peak.o must be of format matrix or sparseMatrix.')
    }

    if (! is(X, 'matrix') & ! is(X, 'sparseMatrix')) {
        stop('X must be of format matrix or sparseMatrix.')
    }

    if (! is(D, 'matrix') & ! is(D, 'sparseMatrix')) {
        stop('D must be of format matrix or sparseMatrix.')
    }

    if (! is(peak.o, 'sparseMatrix')) {
        peak.o <- Matrix(peak.o, sparse = T)
    }

    if (! is(X, 'sparseMatrix')) {
        X <- Matrix(X, sparse = T)
    }

    if (! is(D, 'sparseMatrix')) {
        D <- Matrix(D, sparse = T)
    }

    if (! is(k, 'numeric') & ! is(k, 'integer')) {
        stop('k must be an numeric (integer).')
    }

    if (! is(alpha, 'numeric')) {
        stop('alpha must be of format numeric.')
    }

    if (! is(beta.max.scale, 'numeric')) {
        stop('beta.max.scale must be of format numeric.')
    }

    k <- as.integer(k)

    if (! missing(beta.min)) {
        if (! is(beta.min, numeric)) {
            warning('User must supply a valid beta.min of type numeric, or supply nothing at all. Using default value.')
            beta <- 10.0^seq(beta.max.scale, -3, -1)
        } else {
            beta <- beta.min * 10.0^seq(beta.max.scale, 0, -1)
        }
    } else {
        beta <- 10.0^seq(beta.max.scale, -3, -1)
    }

    # convert to matrix to use NNLM::nnmf function
    message("Running non-negative matrix factorization on X.")
    X.nmf <- nnmf(as.matrix(X), k=k)
    w2 <- X.nmf$W
    h2 <- X.nmf$H

    message("Running non-negative matrix factorization on peak.o.")
    PO.nmf <- nnmf(as.matrix(peak.o), k=k)
    w1 <- PO.nmf$W
    h1 <- PO.nmf$H

    PO.dim <- dim(peak.o)

    mat.init <- initializeMatrix(PO.dim[1],
                          PO.dim[2],
                          dim(X)[2],
                          k,
                          D)

    for (j in 1:length(beta)) {
        dp.ret <- compute_lambda(peak.o,
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
        nmf.ret <- nmf_cluster_sep2_lap(peak.o,
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
    cluster.output <- cluster(H1, H2)

    cluster.output$c1 <- data.frame('ATAC-cluster' = cluster.output$c2)
    cluster.output$c2 <- data.frame('RNA-cluster' = cluster.output$c2)
    if (! missing(rna.cell.label)) {
        if (length(rna.cell.label) != ncol(X)) {
            warning("Cell labels and gene expression matrix are incompatible. Not using cell labels for RNA-seq clustering output.")
        } else {
            rownames(cluster.output$c2) <- rna.cell.label
        }
    }

    if (! missing(atac.cell.label)) {
        if (length(atac.cell.label) != ncol(peak.o)) {
            warning("Cell labels and chromatin accessibility matrix are incompatible. Not using cell labels for ATAC-seq clustering output.")
        } else {
            rownames(cluster.output$c1) <- atac.cell.label
        }
    }

    return(list("W1" = W1,
                "W2" = W2,
                "H1" = H1,
                "H2" = H2,
                "lambda1" = lambda1,
                "lambda2" = lambda2,
                "score" = score,
                "atac_cluster" = cluster.output$c1,
                "rna_cluster" = cluster.output$c2))
}
