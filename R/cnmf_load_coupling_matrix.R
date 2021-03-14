cnmf_load_coupling_matrix <- function(bed.dir,
                                      peak.name,
                                      symbol,
                                      sep.char = '\t',
                                      d0 = 30000,
                                      ...) {
    if (! is(bed.dir, 'character')) {
        stop('bed.dir must be a character (a valid path of the peak_gene_100k_corr.bed file).')
    }

    if (! is.file(bed.dir)) {
        stop('bed.dir is not a valid path to the file.')
    }

    if (! is.character(sep.char) | length(sep.char) > 1) {
        stop('sep.char must be a single character indicating the separater character in the motif target file.')
    }

    if (is (peak.name, "list")) {
        peak.name <- unlist(peak.name, use.names = F)
    }

    if (is (symbol, "list")) {
        symbol <- unlist(symbol, use.names = F)
    }

    if (! is (peak.name, "character")) {
        stop("peak.name must be a vector of character.")
    }

    if (! is (symbol, "character")) {
        stop("symbol must be a vector of character.")
    }

    couple.file <- loadBedFile(bed.dir, sep.char)
    f <- match(couple.file$C1, peak.name, nomatch=0)
    d <- f > 0
    f1 <- match(couple.file$C2, symbol, nomatch=0)
    d1 <- f1 > 0
    couple.file$C4[couple.file$C4 < 0] <- 0
    idx <- which((d * d1) == 1)
    c.array <- exp(-1 / d0 * couple.file$C3[idx]) * couple.file$C4[idx]
    c.array[c.array < 0] <- 0
    D <- sparseMatrix(i=f1[idx],
                      j=f[idx],
                      x=c.array,
                      dims = c(length(symbol), length(peak.name)))
    return(D)
}
