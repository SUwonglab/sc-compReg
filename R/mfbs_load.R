mfbs_load <- function(motif.target.dir,
                      sep.char = '\t',
                      ...) {
    if (! is(motif.target.dir, 'character')) {
        stop('motif.target.dir must be a character (a valid path of the motif target file).')
    }

    if (! is.file(motif.target.dir)) {
        stop('motif.target.dir is not a valid path to the file.')
    }

    if (! is.character(sep.char) | length(sep.char) > 1) {
        stop('sep.char must be a single character indicating the separater character in the motif target file.')
    }

    motif.file <- mfbs.load(motif.target.dir, sep.char)
    return(motif.file)
}
