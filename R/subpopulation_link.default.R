subpopulation.link.default <- function(E.mean.healthy,
                                       E.mean.cll,
                                       O.mean.healthy,
                                       O.mean.cll) {
    if (! is(E.mean.healthy, 'matrix') & ! is(E.mean.healthy, 'Matrix')) {
        stop('E.mean.healthy must be a matrix. Please check your format.')
    }

    if (! is(E.mean.cll, 'matrix') & ! is(E.mean.cll, 'Matrix')) {
        stop('E.mean.cll must be a matrix. Please check your format.')
    }

    if (! is(O.mean.healthy, 'matrix') & ! is(O.mean.healthy, 'Matrix')) {
        stop('O.mean.healthy must be a matrix. Please check your format.')
    }

    if (! is(O.mean.cll, 'matrix') & ! is(O.mean.cll, 'Matrix')) {
        stop('O.mean.cll must be a matrix. Please check your format.')
    }

    output <- subpopulationLink(E.mean.healthy,
                                E.mean.cll,
                                O.mean.healthy,
                                O.mean.cll)

    output$call = this.call

    return(output)
}
