normalize.index <- function(X) {
    return(X - min(X))
}

mfbs.load <- function(motif.path) {
    motif.files <- mfbsLoad(motif.path)
    return(motif.files)
}

maxk <- function(A, k, d) {
    a <- apply(A, d, function(x) {sort(x, decreasing = T)})
    a <- a[, 1:k]
    return(a)
}
