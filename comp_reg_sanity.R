library(scCompReg)
library(R.matlab)
library(Matrix)
library(tictoc)

# rm(list = ls())
s2 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample2.mat')
s1 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample1.mat')

O1=s1$O1
O1.idx = s1$O1.idx
symbol1 = s1$Symbol1
O2 = s2$O2
O2.idx = s2$O2.idx
symbol2 = s2$Symbol2
peak.name1 = s1$PeakName1
E1.idx = s1$E1.idx
E1 = s1$E1
E2 = s2$E2
E2.idx = s2$E2.idx
peak.name2 = s2$PeakName2
peak.name.intersect.dir='/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt'
motif.mat.dir = '/Users/Sophia/Desktop/BioStats/compreg/MotifMatch_human_rmdup.mat'
motif.target.dir ='/Users/Sophia/Desktop/BioStats/compreg/MotifTarget.txt'
peak.gene.prior.dir = '/Users/Sophia/Desktop/BioStats/compreg/peak_gene_prior_intersect.bed'
sep.char=' '

motif.mat <- readMat(motif.mat.dir)
motif.name <- unlist(motif.mat$motifName, use.names=F)
motif.weight <- as.numeric(unlist(motif.mat$motifWeight, use.names=F))
match2 <- unlist(motif.mat$Match2, use.names = F)
m2.half.idx <- length(match2) / 2
match2 <- list('a' = match2[1: m2.half.idx],
               'b' = match2[(m2.half.idx + 1):length(match2)])


output = sc_compreg(O1,
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
           '/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt',
           '/Users/Sophia/Desktop/BioStats/compreg/MotifTarget.txt',
           '/Users/Sophia/Desktop/BioStats/compreg/peak_gene_prior_intersect.bed')
write.table(output$hub.tf[[2]], 'tf2.txt')
write.table(output$diff.net[[2]], 'diff_net2.txt')

pni = comp_reg_preprocess('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', token='\t')

symbol =clust.profile.output$symbol
symbol = unlist(symbol, use.names=F)
tf.binding = mfbs.output$tf.binding
tf.name = mfbs.output$tf.name
elem.name =clust.profile.output$elem.name
peak.gene.prior.path =peak.gene.prior.dir
E1 = s1$E1
E1.idx=s1$E1.idx
E2 =s2$E2
E2.idx =s2$E2.idx
O1.mean=clust.profile.output$O1.mean
O2.mean=clust.profile.output$O2.mean
match.mat=subpop.link.output$match
thresh = 0.2
sig.level = 0.05
num.top.tf = 5000
d0.default = 500000
verbose=T


# tic()
output = cluster.profile(Matrix(s1$O1, sparse=T),
               Matrix(s1$E1, sparse=T),
                as.integer(s1$O1.idx),
               as.integer(s1$E1.idx),
               as.vector(unlist(s1$Symbol1, use.names = F)),
                as.vector(unlist(s1$PeakName1, use.names = F)),
               Matrix(s2$O2, sparse=T),
               Matrix(s2$E2, sparse=T),
               as.integer(s2$O2.idx),
               as.integer(s2$E2.idx),
               as.vector(unlist(s2$Symbol2, use.names = F)),
               as.vector(unlist(s2$PeakName2, use.names = F)),
                pni$vo,
               pni$vt
               )
# toc()

O1 = s1$O1
E1 = s1$E1
E1.idx = s1$E1.idx
O1.idx = s1$O1.idx
symb1 = s1$Symbol1
symb2 = s2$Symbol2
O2 = s2$O2
E2 = s2$E2
E2.idx = s2$E2.idx
O2.idx = s2$O2.idx
pk.name1 = unlist(s1$PeakName1, use.names = F)
pk.name2 = unlist(s2$PeakName2, use.names = F)
pk.name.intersect1 = pni$vo
pk.name.intersect2 = pni$vt

K1 = max(max(O1.idx), max(E1.idx))
symbol = intersect(symb1, symb2)
f1 = match(symbol, symb1)
E1.mean = matrix(NA, length(f1), K1)
for (i in 1:K1) {
    gp = which(E1.idx == i)
    if (length(gp) == 0) next
    E1.mean[, i] = rowMeans(E1[f1, gp])
}



K2 = max(max(O2.idx), max(E2.idx))
f2 = match(symbol, symb2)
E2.mean = matrix(NA, length(f2), K2)
for (i in 1:K2) {
    gp = which(E2.idx == i)
    if (length(gp) == 0) next
    E2.mean[, i] = rowMeans(E2[f2, gp])
}

pk.name1 = unlist(pk.name1, use.names = F)
pk.name2 = unlist(pk.name2, use.names=F)

pf1 = match(pk.name.intersect1, pk.name1, nomatch=0)
pf2 = match(pk.name.intersect2, pk.name2, nomatch=0)
d1 = is.element(pk.name1, pk.name.intersect1)
d2 = is.element(pk.name2, pk.name.intersect2)
print(length(d1))
print(length(d2))

elem.name = c(pk.name.intersect1, pk.name1[d1 == 0],
                  pk.name2[d2 == 0])

m = length(pk.name.intersect1)
elem.len = length(elem.name)
O1.mean = matrix(0, elem.len, K1)
O2.mean = matrix(0, elem.len, K2)

sum.d1 <- sum(d1 == 0)

for (i in 1:K1) {
    gp = which(O1.idx == i)
    if (length(gp) == 0) next
    O1.mean[1:m, i] = rowMeans(O1[pf1, gp] > 0)
    O1.mean[(1+m) : (m + sum.d1), i] = rowMeans(O1[d1 == 0, gp] > 0)
}

for (i in 1:K2) {
    gp = which(O2.idx == i)
    if (length(gp) == 0) next
    O2.mean[1:m, i] = rowMeans(O2[pf2, gp] > 0)
    O2.mean[(m + sum.d1 + 1) : elem.len, i] = rowMeans(O2[d2 == 0, gp] > 0)
}


O1.mean = sweep(O1.mean, 2, colMeans(O1.mean), '/')
O2.mean = sweep(O2.mean, 2, colMeans(O2.mean), '/')

K1 = max(max(s1$O1.idx), max(s1$E1.idx))
K2 = max(max(s2$O2.idx), max(s2$E2.idx))
symbol = intersect(s1$Symbol1, s2$Symbol2)
f1 = match(symbol, s1$Symbol1)
f2 = match(symbol, s2$Symbol2)
E1.mean = matrix(NA, length(f1), K1)
for (i in 1:K1) {
    E1.mean[, i] = apply(s1$E1[f1, s1$E1.idx == i], 1, mean)
}



new.output = subpopulation.link(output$E1.mean,
                                output$E2.mean,
                                output$O1.mean,
                                output$O2.mean)
# write.csv(output$O1.mean, 'O1mean.csv')
# write.csv(output$O2.mean, 'O2mean.csv')
# write.csv(output$E1.mean, 'E1mean.csv')
# write.csv(output$E2.mean, 'E2mean.csv')

pkname = read.table('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', header=F)
length(pkname$V1)
length(pkname$V2)

E.mean.healthy = E1.mean
E.mean.cll = E2.mean
O.mean.healthy = O1.mean
O.mean.cll = O2.mean

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



m = readMat('/Users/Sophia/Desktop/BioStats/compreg/MotifMatch_human_rmdup.mat')




# tic()
motif.file <- mfbs.load('/Users/Sophia/Desktop/BioStats/compreg/MotifTarget.txt')
# toc()

f3 <- motif.file$C3
f1 <- match(motif.file$C1, elem.name, nomatch = 0)
d1 <- f1 > 0

motif.name <- unlist(m$motifName, use.names=F)
motif.weight <- as.numeric(unlist(m$motifWeight, use.names=F))

f2 <- match(motif.file$C2, motif.name, nomatch = 0)
d2 <- f2 > 0

t1 <- setdiff(seq(1, length(motif.name), 1), unique(f2))


f2 <- c(f2[(d1 * d2) == 1], t1)
t1.len <- length(t1)
f1 <- c(f1[(d1 * d2) == 1], rep(1, t1.len))
f3 <- c(f3[(d1 * d2) == 1], rep(0, t1.len))

t1 <- setdiff(seq(1, length(elem.name), 1), unique(f1))
f1 <- c(f1, t1)
t1.len <- length(t1)
f2 <- c(f2, rep(1, t1.len))
f3 <- c(f3, rep(0, t1.len))


symbol <- unlist(symbol, use.names=F)

motif.binding <- sparseMatrix(dims = c(length(motif.name),length(elem.name)),
                              i = f2,
                              j = f1,
                              x = f3)


motif.binding <- mult(Matrix(diag(1 / (motif.weight + 0.1)), sparse=T), motif.binding)
motif.binding@x <- log(1 + motif.binding@x)

match2 <- unlist(m$Match2, use.names = F)

m2.half.idx <- length(match2) / 2
match2 <- list('a' = match2[1: m2.half.idx],
               'b' = match2[(m2.half.idx + 1):length(match2)])

tf.name <- intersect(symbol, unique(match2$b))
tf.name.len <- length(tf.name)
tf.binding <- matrix(0, tf.name.len, length(elem.name))

mf1 <- match(match2$a, motif.name, nomatch=0)
mf2 <- match(match2$b, tf.name, nomatch=0)

for (i in 1:tf.name.len) {
    a <- which(mf2 == i)
    if (length(a) > 1) {
        tf.binding[i, ] <- colMax(motif.binding[mf1[a], ])
    } else if (length(a) == 1) {
        tf.binding[i, ] <- motif.binding[mf1[a], ]
    }
}

# tf.binding <- Matrix(tf.binding, sparse = T)



############# comp reg #############
symbol = unlist(symbol, use.names = F)
a <- threshK(tf.binding, 5000)


tf.binding[tf.binding - a[, ncol(a)] < 0] <- 0

file <- compRegLoad('/Users/Sophia/Desktop/BioStats/compreg/peak_gene_prior_intersect.bed')
f <- match(file$C1, elem.name, nomatch=0)
d <- f > 0
f1 <- match(file$C2, symbol, nomatch=0)
d1 <- f1 > 0
indices <- (d * d1) == 1
fc <- cbind(f[indices], f1[indices])

f2 <- unique(fc, orient = 'r')
ic <- match(do.call(paste, data.frame(fc)), do.call(paste, data.frame(f2)), nomatch=0)
ic <- as.numeric((ic - 1))

c3 <- accumArrayMin(ic, as.numeric(file$C3[indices]))
c4 <- accumArrayMin(ic, as.numeric(file$C4[indices]))


thresh <- 0.2
c4[c4 < thresh] <- 0
d0 <- 500000
c <- as.numeric(exp(-1 * c3 / d0) * c4)
beta.sparse <- sparseMatrix(dims = c(length(elem.name), length(symbol)),
                     j = f2[, 2],
                     i = f2[, 1],
                     x = c)

diff.net <- list()
hub.tf <- list()




tic()
a = sweep(tf.binding, 2, O1.mean[, i1], '*') %*% beta.sparse
toc()
sp = sweep(tf.binding, 2, O1.mean[, i1], '*')

tic()
a = computeBZ(tf.binding, O1.mean[, i1], beta.sparse)
toc()


for (ii in 1:nrow(match.mat)) {
    i1 <- match.mat[ii, 1]
    i2 <- match.mat[ii, 2]
    B01 <- sweep(tf.binding, 2, O1.mean[, i1], '*') %*% beta.sparse
    B02 <- sweep(tf.binding, 2, O2.mean[, i2], '*') %*% beta.sparse
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

    # adjusted.p.val <- read.table('/Users/Sophia/Desktop/BioStats/scCompReg/adj_p.txt',
    #                              header=F,
    #                              sep='\t')
    diff.gene <- which(adjusted.p.val < p.val.thresh)



    corr.test.stat1 <- corrTest(t(TF[, 1:n1]), TG1)
    p1 <- 2 * pt(abs(corr.test.stat1), df = nrow(TG1) - 2, lower.tail = F)
    corr.test.stat2 <- corrTest(t(TF[, (1 + n1) : (n1 + n2)]), TG2)
    p2 <- 2 * pt(abs(corr.test.stat2), df = nrow(TG2) - 2, lower.tail = F)
    p.combine <- pmin(p1, p2, na.rm = T)
    net.idx <- which(p.combine < 0.05, arr.ind = T)
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
            LR.summary.id <- c(LR.summary.id, LR.summary.id[i])
            LR.summary.diff.gene <- c(LR.summary.diff.gene, diff.gene.j)
            LR.summary.stat <- c(LR.summary.stat, lr.ret$stat)
            LR.summary.p.val <- c(LR.summary.p.val, lr.ret$p)
        }
    }
    LR.summary <- cbind(LR.summary.id,
                        LR.summary.diff.gene,
                        LR.summary.stat,
                        LR.summary.p.val)

    LR.summary <- LR.summary[which(!is.na(LR.summary.stat)
                                   & !is.infinite(LR.summary.stat)), ]

    LR.summary[LR.summary[, 3] < 10^(-16), 3] <- 10^(-16)

    p.hat <- fit.gamma.quantile.matching(LR.summary[, 3],
                                         seq(0.01, 0.2, 0.01))
    p <- 1 - pgamma(LR.summary[, 3], p.hat[1], p.hat[2])
    adj.p <- p.adjust(p, "fdr")

    LR.summary <- cbind(LR.summary,
                        p,
                        adj.p)
    idx <- which(adj.p < 0.1)
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







