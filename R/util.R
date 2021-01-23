normalize.index <- function(X) {
    return(X - min(X))
}

mfbs.load <- function(motif.path) {
    motif.files <- mfbsLoad(motif.path)
    return(motif.files)
}
