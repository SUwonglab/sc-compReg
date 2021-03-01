mfbs_load <- function(motif.target.dir) {
    if (! is(motif.target.dir, 'character')) {
        stop('motif.target.dir must be a character (a valid path of the motif target file).')
    }

    if (! is.file(motif.target.dir)) {
        stop('motif.target.dir is not a valid path to the file.')
    }

    motif.file <- mfbs.load(motif.target.dir)
    return(motif.file)
}
