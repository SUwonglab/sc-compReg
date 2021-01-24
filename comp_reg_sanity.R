library(scCompReg)
library(R.matlab)
library(Matrix)
library(tictoc)

pni = comp_reg_preprocess('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', token='\t')
s2 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample2.mat')
s1 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample1.mat')

tic()
output = cluster_profile(Matrix(s1$O1, sparse=T),
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
toc()

O1 = s1$O1
E1 = s1$E1
E1.idx = s1$E1.idx
O1.idx = s1$O1.idx
symb1 = s1$Symbol1
symb2 = s2$Symbol2
O2 = s2$O2
E2 = s2$E2
E2.idx = s2$E2.idx
O2.idx = s2$E2.idx
pk.name1 = s1$PeakName1
pk.name2 = s2$PeakName2
pk.name.intersect1 = pni$vo
pk.name.intersect2 = pni$vt

K1 = max(max(O1.idx), max(E1.idx))
symbol = intersect(symb1, symb2)
f1 = match(symbol, symb1)
E1.mean = matrix(NA, length(f1), K1)
for (i in 1:K1) {
    gp = which(E1.idx == i)
    if (length(gp) == 0) next
    E1.mean[, i] = Matrix::rowMeans(E1[f1, gp])
}



K2 = max(max(O2.idx), max(E2.idx))
f2 = match(symbol, symb2)
E2.mean = matrix(NA, length(f2), K2)
for (i in 1:K2) {
    gp = which(E2.idx == i)
    if (length(gp) == 0) next
    E2.mean[, i] = Matrix::rowMeans(E2[f2, gp])
}

pf1 = match(pk.name.intersect1, pk.name1)
pf2 = match(pk.name.intersect2, pk.name2)
d1 = is.element(pk.name1, pk.name.intersect1)
d2 = is.element(pk.name2, pk.name.intersect2)
print(length(d1))
print(length(d2))
pk.name1 = unlist(pk.name1, use.names = F)
pk.name2 = unlist(pk.name2, use.names=F)
elem.name = c(pk.name.intersect1, pk.name1[d1 == 0],
                  pk.name2[d2 == 0])

m = length(pk.name.intersect1)
elem.len = length(elem.name)
O1.mean = matrix(NA, elem.len, K1)
O2.mean = matrix(NA, elem.len, K2)

for (i in 1:K1) {
    gp = which(O1.idx == i)
    if (length(gp) == 0) next
    O1.mean[1:m, i] = Matrix::rowMeans(O1[pf1, gp] > 0)
    O1.mean[(1+m) : (m + sum(d1 == 0)), i] = Matrix::rowMeans(O1[d1 == 0, gp] > 0)
}

for (i in 1:K2) {
    gp = which(O2.idx == i)
    if (length(gp) == 0) next
    O2.mean[1:m, i] = Matrix::rowMeans(O2[pf2, gp] > 0)
    O2.mean[(m + sum(d1 == 0) + 1) : elem.len, i] = Matrix::rowMeans(O2[d2 == 0, gp] > 0)
}

O1.mean = O1.mean / Matrix::colMeans(O1.mean)
O2.mean = O2.mean / Matrix::colMeans(O2.mean)


K1 = max(max(s1$O1.idx), max(s1$E1.idx))
K2 = max(max(s2$O2.idx), max(s2$E2.idx))
symbol = intersect(s1$Symbol1, s2$Symbol2)
f1 = match(symbol, s1$Symbol1)
f2 = match(symbol, s2$Symbol2)
E1.mean = matrix(NA, length(f1), K1)
for (i in 1:K1) {
    E1.mean[, i] = apply(s1$E1[f1, s1$E1.idx == i], 1, mean)
}


output$E1.mean == E1.mean

new.output = subpopulation.link(output$E1.mean,
                                output$E2.mean,
                                output$O1.mean,
                                output$O2.mean)
write.csv(output$O1.mean, 'O1mean.csv')
write.csv(output$O2.mean, 'O2mean.csv')
write.csv(output$E1.mean, 'E1mean.csv')
write.csv(output$E2.mean, 'E2mean.csv')

pkname = read.table('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', header=F)
length(pkname$V1)
length(pkname$V2)

E.mean.healthy = output$E1.mean
E.mean.cll = output$E2.mean
O.mean.healthy = output$O1.mean
O.mean.cll = output$O2.mean

r1 <- cor(E.mean.healthy, E.mean.cll)
r2 <- cor(O.mean.healthy, O.mean.cll)
rr1 <- r1 - colSums(r1) %*% t(rowSums(r1)) / sum(r1)
rr2 <- r2 - colSums(r2) %*% t(rowSums(r2)) / sum(r2)
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
match <- matrix(NA, a.dim, 7)
S1 <- matrix(NA, a.dim, 1)
S2 <- matrix(NA, a.dim, 1)
match.idx <- 1
for (i in 1:a.dim) {
    idx <- is.element(a[i, 1], S1) + is.element(a[i, 2], S2)
    S1[i] <- a[i, 1]
    S2[i] <- a[i, 2]
    if (idx < 2) {
        match[match.idx, ] <- a[i, ]
        match.idx = match.idx + 1
    }
}

match <- match[rowSums(is.na(match)) != ncol(match), ]
write.csv(match, 'match.csv')

output = list()
output$match <- match
output$call <- this.call






m = readMat('/Users/Sophia/Desktop/BioStats/compreg/MotifMatch_human_rmdup.mat')




tic()
motif.file <- mfbs.load('/Users/Sophia/Desktop/BioStats/compreg/MotifTarget.txt')
toc()

f3 <- motif.file$C3
d1 <- is.element(motif.file$C1, elem.name)
f1 <- match(motif.file$C1, elem.name, nomatch = 0)
motif.name <- unlist(m$motifName, use.names=F)
motif.weight <- as.numeric(unlist(m$motifWeight, use.names=F))

d2 <- is.element(motif.file$C2, motif.name)
f2 <- match(motif.file$C2, motif.name, nomatch = 0)

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


motif.binding <- sparseMatrix(dims = c(length(motif.name),length(elem.name)),
                              i = f2,
                              j = f1,
                              x = f3)
motif.weight.len <- length(motif.weight)


motif.binding <- mult(Matrix(diag(1 / (motif.weight + 0.1),
                      nrow = motif.weight.len,
                      ncol = motif.weight.len), sparse=T), motif.binding)
motif.binding@x <- log(motif.binding@x)

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
    print(i)
    a <- which(mf2 == i)
    if (length(a) > 1) {
        tf.binding[i, ] <- max(motif.binding[mf1[a], ])
    } else if (length(a) == 1) {
        tf.binding[i, ] <- motif.binding[mf1[a], ]
    }
}

tf.binding <- Matrix(tf.binding, sparse = T)


