sc_compreg.default <- function(O1,
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
                               motif.name,
                               motif.weight,
                               match2,
                               peak.name.intersect.dir,
                               motif.mat.dir,
                               motif.target.dir,
                               peak.gene.prior.dir,
                               sep.char = '\t',
                               thresh = 0.2,
                               sig.level = 0.05,
                               num.top.tf = 5000,
                               d0.default = 500000,
                               verbose = TRUE) {
    if (! is(O1, 'sparseMatrix')) {
        if (! is(O1, 'matrix')) {
            stop('O1 must be a matrix or sparseMatrix. Please check your format.')
        }
        O1 <- Matrix(O1, sparse = T)
    }

    if (! is(O2, 'sparseMatrix')) {
        if (! is(O2, 'matrix')) {
            stop('O2 must be a matrix or sparseMatrix. Please check your format.')
        }
        O2 <- Matrix(O2, sparse = T)
    }

    if (! is(E1, 'sparseMatrix')) {
        if (! is(E1, 'matrix')) {
            stop('E1 must be a matrix or sparseMatrix. Please check your format.')
        }
        E1 <- Matrix(E1, sparse = T)
    }

    if (! is(E2, 'sparseMatrix')) {
        if ( ! is(E2, 'matrix')) {
            stop('E2 must be a matrix or sparseMatrix. Please check your format.')
        }
        E2 <- Matrix(E2, sparse = T)
    }

    if (! is(O1.idx, 'numeric') &
        ! is(O1.idx, 'matrix') &
        ! is(O1.idx, 'vector') &
        ! is(E1.idx, 'integer')) {
        if (length(O1.idx < 1) | O1.idx[1] %% 1 != 0) {
            stop('O1.idx must be a vector of integer (indices).')
        }
    }

    O1.idx <- as.integer(O1.idx)

    if (! is(E1.idx, 'numeric') &
        ! is(E1.idx, 'matrix') &
        ! is(E1.idx, 'vector') &
        ! is(E1.idx, 'integer')) {
        if (length(E1.idx < 1) | E1.idx[1] %% 1 != 0) {
            stop('E1.idx must be a vector of integer (indices).')
        }
    }

    E1.idx <- as.integer(E1.idx)

    if (! is(symbol1, 'vector')) {
        stop('symbol1 must be a vector of characters.')
    }

    if (is(peak.name1, 'list')) {
        peak.name1 <- unlist(peak.name1, use.names = F)
    }

    if (! is(peak.name1, 'vector')) {
        stop('peak.name1 must be a vector of characters.')
    }


    if (! is(O2.idx, 'numeric') &
        ! is(O2.idx, 'matrix') &
        ! is(O2.idx, 'vector') &
        ! is(E1.idx, 'integer')) {
        if (length(O2.idx < 1) | O2.idx[1] %% 1 != 0) {
            stop('O2.idx must be a vector of integer (indices).')
        }
    }

    O2.idx <- as.integer(O2.idx)

    if (! is(E2.idx, 'numeric') &
        ! is(E2.idx, 'matrix') &
        ! is(E2.idx, 'vector') &
        ! is(E1.idx, 'integer')) {
        if (length(E2.idx < 1) | E2.idx[1] %% 1 != 0) {
            stop('E2.idx must be a vector of integer (indices).')
        }
    }

    E2.idx <- as.integer(E2.idx)

    if (! is(symbol2, 'vector')) {
        stop('symbol2 must be a vector of characters.')
    }

    if (is(peak.name2, 'list')) {
        peak.name2 <- unlist(peak.name2, use.names = F)
    }

    if (! is(peak.name2, 'vector')) {
        stop('peak.name2 must be a vector of characters.')
    }

    pni.file <- comp_reg_preprocess(peak.name.intersect.dir, token=sep.char)

    clust.profile.output <- cluster_profile(O1,
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
                                            pni.file$vo,
                                            pni.file$vt)

    subpop.link.output <- subpopulation_link(clust.profile.output$E1.mean,
                                             clust.profile.output$E2.mean,
                                             clust.profile.output$O1.mean,
                                             clust.profile.output$O2.mean)

    clust.profile.output$symbol <- unlist(clust.profile.output$symbol, use.names = F)

    motif.file <- mfbs.load(motif.target.dir)
    mfbs.output <- mfbs(clust.profile.output$elem.name,
                        clust.profile.output$symbol,
                        motif.name,
                        motif.weight,
                        match2,
                        motif.file)

    compreg.output <- compreg(clust.profile.output$symbol,
                              mfbs.output$tf.binding,
                              mfbs.output$tf.name,
                              clust.profile.output$elem.name,
                              peak.gene.prior.dir,
                              s1$E1,
                              s1$E1.idx,
                              s2$E2,
                              s2$E2.idx,
                              clust.profile.output$O1.mean,
                              clust.profile.output$O2.mean,
                              subpop.link.output$match,
                              thresh,
                              sig.level,
                              num.top.tf,
                              d0.default,
                              verbose)

    compreg.output$call <- match.call()
    return(compreg.output)
}
