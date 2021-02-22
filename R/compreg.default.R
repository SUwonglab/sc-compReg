compreg.default <- function(symbol,
                            tf.binding,
                            tf.name,
                            elem.name,
                            peak.gene.prior.path,
                            E1,
                            E1.idx,
                            E2,
                            E2.idx,
                            O1.mean,
                            O2.mean,
                            match.mat,
                            thresh,
                            sig.level,
                            num.top.tf,
                            d0.default,
                            verbose) {
    if (!is(peak.gene.prior.path, 'character')) {
        stop('peak.gene.prior.path must be of class character.')
    }

    this.call <- match.call()

    a <- threshK(tf.binding, num.top.tf)

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

    # adjust for index starting off at one in R to call function in cpp
    ic <- as.numeric((ic - 1))

    c3 <- accumArrayMin(ic, as.numeric(file$C3[indices]))
    c4 <- accumArrayMin(ic, as.numeric(file$C4[indices]))

    c4[c4 < thresh] <- 0
    c <- as.numeric(exp(-1 * c3 / d0.default) * c4)
    beta.sparse <- sparseMatrix(dims = c(length(elem.name), length(symbol)),
                                j = f2[, 2],
                                i = f2[, 1],
                                x = c)


    diff.net <- list()
    hub.tf <- list()

    for (ii in 1:nrow(match.mat)) {
        if (verbose) {
            message(paste('Running', ii, 'th link subpopulation.'))
        }
        i1 <- match.mat[ii, 1]
        i2 <- match.mat[ii, 2]
        B01 <- computeBZ(tf.binding, O1.mean[, i1], beta.sparse)
        B02 <- computeBZ(tf.binding, O2.mean[, i2], beta.sparse)
        TG1 <- E1[, E1.idx == i1]
        TG2 <- E2[, E2.idx == i2]
        TG2 <- TG2 * mean(TG1) / mean(TG2)
        f <- match(tf.name, symbol, nomatch=0)
        TF <- cbind(TG1[f, ], TG2[f, ])
        n1 <- ncol(TG1)
        n2 <- ncol(TG2)
        TG1 <- t(TG1)
        TG2 <- t(TG2)
        p.val <- sapply(1:ncol(TG1),
                        function(x) t.test(TG1[,x], TG2[,x], var.equal = T)$p.value)
        adjusted.p.val <- p.adjust(p.val, "fdr")

        diff.gene <- which(adjusted.p.val < sig.level)

        corr.test.stat1 <- corrTest(t(TF[, 1:n1]), TG1)
        p1 <- 2 * pt(abs(corr.test.stat1), df = nrow(TG1) - 2, lower.tail = F)
        corr.test.stat2 <- corrTest(t(TF[, (1 + n1) : (n1 + n2)]), TG2)
        p2 <- 2 * pt(abs(corr.test.stat2), df = nrow(TG2) - 2, lower.tail = F)
        p.combine <- pmin(p1, p2, na.rm = T)
        net.idx <- which(p.combine < sig.level, arr.ind = T)
        LR.summary.id <- c()
        LR.summary.diff.gene <- c()
        LR.summary.stat <- c()
        LR.summary.p.val <- c()

        for (j in 1:length(diff.gene)) {
            OTF1 <- t(B01[, diff.gene[j]] * TF[, 1:n1])
            OTF2 <- t(B02[, diff.gene[j]] * TF[, (1 + n1) : (n1 + n2)])
            diff.gene.j <- diff.gene[j]
            id1 <- net.idx[net.idx[, 2] == diff.gene.j, 1]
            id <- which((colSums(abs(OTF1)) + colSums(abs(OTF2))) > 0)
            id <- intersect(id, id1)
            if (length(id) == 0) next

            for (i in 1:length(id)) {
                X1 <- cbind(OTF1[, id[i]], TG1[, diff.gene.j])
                X2 <- cbind(OTF2[, id[i]], TG2[, diff.gene.j])
                lr.ret <- bivariate.normal.conditional.lr(X1, X2)
                LR.summary.id <- c(LR.summary.id, id[i])
                LR.summary.diff.gene <- c(LR.summary.diff.gene, diff.gene.j)
                LR.summary.stat <- c(LR.summary.stat, lr.ret$stat)
                LR.summary.p.val <- c(LR.summary.p.val, lr.ret$p)
            }
        }
        LR.summary <- cbind(LR.summary.id,
                            LR.summary.diff.gene,
                            LR.summary.stat,
                            LR.summary.p.val)

        LR.summary <- LR.summary[which((!is.na(LR.summary.stat)) & (!is.infinite(LR.summary.stat))), ]

        LR.summary[LR.summary[, 3] < 10^(-16), 3] <- 10^(-16)

        p.hat <- fit.gamma.quantile.matching(LR.summary[, 3],
                                             seq(0.01, 0.2, 0.01))
        p <- 1 - pgamma(LR.summary[, 3], p.hat[1], p.hat[2])
        adj.p <- p.adjust(p, "fdr")

        LR.summary <- cbind(LR.summary,
                            p,
                            adj.p)
        idx <- which(adj.p < (sig.level*2))
        sorted.idx <- order(LR.summary[idx, 3], decreasing = T)
        idx1 <- idx[f]

        diff.net.ii <- cbind(tf.name[LR.summary[idx1, 1]],
                             symbol[LR.summary[idx1, 2]],
                             LR.summary[idx1, 3:6])
        u <- unique(diff.net.ii[, 1])
        u.f <- match(diff.net.ii[, 1], u, nomatch = 0)
        diff.net[[ii]] <- diff.net.ii
        u.len <- length(u)
        u.count <- numeric(u.len)
        for (i in 1:u.len) {
            u.count[i] <- sum(f == i)
        }
        u.f1 <- order(u.count, decreasing=T)
        u.d1 <- u.count[u.f1]
        hub.tf[[ii]] <- cbind(u[u.f1], u.d1)
    }

    output <- list()
    output$call <- this.call
    output$hub.tf <- hub.tf
    output$diff.net <- diff.net
    return(output)
}
