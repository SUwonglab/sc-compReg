cluster.profile.default <- function(O1,
                            E1,
                            O1.idx,
                            E1.idx,
                            symb1,
                            pk.name1,
                            O2,
                            E2,
                            O2.idx,
                            E2.idx,
                            symb2,
                            pk.name2,
                            pk.name.intersect1,
                            pk.name.intersect2) {

    if (! is(O1, 'sparseMatrix')) {
        if (! is(O1, 'matrix')) {
            stop('O1 must be a matrix or sparseMatrix. Please check your format.')
        }
        O1 = Matrix(O1, sparse = T)
    }

    if (! is(O2, 'sparseMatrix')) {
        if (! is(O2, 'matrix')) {
            stop('O2 must be a matrix or sparseMatrix. Please check your format.')
        }
        O2 = Matrix(O2, sparse = T)
    }

    if (! is(E1, 'sparseMatrix')) {
        if (! is(E1, 'matrix')) {
            stop('E1 must be a matrix or sparseMatrix. Please check your format.')
        }
        E1 = Matrix(E1, sparse = T)
    }

    if (! is(E2, 'sparseMatrix')) {
        if ( ! is(E2, 'matrix')) {
            stop('E2 must be a matrix or sparseMatrix. Please check your format.')
        }
        E2 = Matrix(E2, sparse = T)
    }

    if (! is(O1.idx, 'numeric')) {
        stop('O1.idx must be a vector of integer (indices).')
    }

    O1.idx = as.integer(O1.idx)
    O1.idx = normalize.index(O1.idx) #scale down with starting index of 0

    if (! is(E1.idx, 'numeric') | ! is(E1.idx, 'vector')) {
        stop('E1.idx must be a vector of integer (indices).')
    }

    E1.idx = as.integer(E1.idx)
    E1.idx = normalize.index(E1.idx) #scale down to have index starting at 0

    if (! is(symb1, 'vector')) {
        stop('symb1 must be a vector of characters.')
    }

    if (! is(pk.name1, 'vector')) {
        stop('pk.name1 must be a vector of characters.')
    }


    if (! is(O2.idx, 'numeric') | ! is(O2.idx, 'vector')) {
        stop('O2.idx must be a vector of integer (indices).')
    }

    O2.idx = as.integer(O2.idx)
    O2.idx = normalize.index(O2.idx) #scale down with starting index of 0

    if (! is(E2.idx, 'numeric') | ! is(E2.idx, 'vector')) {
        stop('E2.idx must be a vector of integer (indices).')
    }

    E2.idx = as.integer(E2.idx)
    E2.idx = normalize.index(E2.idx) #scale down to have index starting at 0

    if (! is(symb2, 'vector')) {
        stop('symb2 must be a vector of characters.')
    }

    if (! is(pk.name2, 'vector')) {
        stop('pk.name2 must be a vector of characters.')
    }

    if (! is(pk.name.intersect1, 'vector')) {
        stop('pk.name.intersect1 must be a vector of characters.')
    }

    if (! is(pk.name.intersect2, 'vector')) {
        stop('pk.name.intersect2 must be a vector of characters.')
    }

    this.call = match.call()

    output = clusterProfile(O1,
                            E1,
                            O1.idx,
                            E1.idx,
                            symb1,
                            pk.name1,
                            O2,
                            E2,
                            O2.idx,
                            E2.idx,
                            symb2,
                            pk.name2,
                            pk.name.intersect1,
                            pk.name.intersect2)
    names(output) = c('E1.mean', 'E2.mean', 'symbol', 'O1.mean', 'O2.mean',
                      'elem.name')
    output$call = this.call

    return(output)
}
