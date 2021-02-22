normalize.index <- function(X) {
    return(X - min(X))
}

mfbs.load <- function(motif.path) {
    motif.files <- mfbsLoad(motif.path)
    return(motif.files)
}

bivariate.normal.conditional.lr <- function(X1, X2,
                                            sigma.min.pseudocount = 0) {
    XX <- rbind(X1, X2)
    beta <- pinv(cbind(rep(1, nrow(XX)), XX[, 2])) %*% XX[, 1]
    est.X <- beta[1] + beta[2] * XX[, 2]
    cond.var <- var(XX[, 1] - est.X)
    L.M0 <- sum(log(dnorm(XX[, 1],
                          mean = est.X,
                          sd = sqrt(cond.var) + sigma.min.pseudocount)))
    beta1 <- pinv(cbind(rep(1, nrow(X1)), X1[, 2])) %*% X1[, 1]
    est.X1 <- beta1[1] + beta1[2] * X1[, 2]
    cond.var1 <- var(X1[, 1] - est.X1)
    beta2 <- pinv(cbind(rep(1, nrow(X2)), X2[, 2])) %*% X2[, 1]
    est.X2 <- beta2[1] + beta2[2] * X2[, 2]
    cond.var2 <- var(X2[, 1] - est.X2)
    L.M1.1 <- sum(log(dnorm(X1[, 1],
                          mean = est.X1,
                          sd = sqrt(cond.var1) + sigma.min.pseudocount)))
    L.M1.2 <- sum(log(dnorm(X2[, 1],
                          mean = est.X2,
                          sd = sqrt(cond.var2) + sigma.min.pseudocount)))
    L.M1 <- L.M1.1 + L.M1.2
    stat <- 2 * (L.M1 - L.M0)
    return(list('p' = 1 - pchisq(stat, df = 3),
                'stat' = stat))
}

fit.gamma.quantile.matching <- function(X,
                                        eta) {
    cut <- quantile(X, probs = eta, names=F)
    fun <- function(alpha) {
        alpha[2] <- abs(alpha[2])
        return(sum(
            (pgamma(cut,
                   shape = alpha[1],
                   rate = (1 / alpha[2])) - eta) ^ 2 ) )
    }

    cut.end <- cut[length(cut)]
    gamma.cut <- sapply(X,
                        function(z) {min(z, cut.end)}
                        )
    x0 <- egamma(gamma.cut)$parameters
    x <- fminunc(x0, fun,
                maxiter = 600)
    x$par <- abs(x$par)
    return(x$par)
}














