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

    K1 = max(max(O1.idx), max(E1.idx))
    symbol = intersect(Symbol1, Symbol2)
    f1 = match(symbol, Symbol1)
    E1.mean = matrix(NA, length(f1), K1)
    for (i in 1:K1) {
        gp = which(E1.idx == i)
        if (length(gp) == 0) next
        E1.mean[, i] = rowMeans(E1[f1, gp])
    }

    K2 = max(max(O2.idx), max(E2.idx))
    f2 = match(symbol, Symbol2)
    E2.mean = matrix(NA, length(f2), K2)
    for (i in 1:K2) {
        gp = which(E2.idx == i)
        if (length(gp) == 0) next
        E2.mean[, i] = rowMeans(E2[f2, gp])
    }

    pf1 = match(pk.name.intersect1, pk.name1)
    pf2 = match(pk.name.intersect2, pk.name2)
    d1 = match(pk.name1, pk.name.intersect1)
    d2 = match(pk.name2, pk.name.intersect2)
    elem.name = rbind(pk.name.intersect1, pk.name1[d1 == 0],
                      pk.name2[d2 == 0])
    m = length(pk.name.intersect1)
    elem.len = length(elem.name)
    O1.mean = matrix(NA, elem.len, K1)
    O2.mean = matrix(NA, elem.len, K2)

    for (i in 1:K1) {
        gp = which(O1.idx == i)
        if (length(gp) == 0) next
        O1.mean[1:m, i] = rowMeans(O1[pf1, gp] > 0)
        O1.mean[(1+m):sum(d1 == 0), i] = rowMeans(O1[d1 == 0, gp] > 0)
    }

    for (i in 1:K2) {
        gp = which(O2.idx == i)
        if (length(gp) == 0) next
        O2.mean[1:m, i] = rowMeans(O2[f2, gp] > 0)
        O2.mean[(m + sum(d1 == 0) + 1) : elem.len, i] = rowMeans(O2[f2 == 0, gp] > 0)
    }

    O1.mean = O1.mean / colMeans(O1.mean)
    O2.mean = O2.mean / colMeans(O2.mean)

    output = list('E1.mean' = E1.mean,
                  'E2.mean' = E2.mean,
                  'symbol' = symbol,
                  'O1.mean' = O1.mean,
                  'O2.mean' = O2.mean,
                  'elem.name' = elem.name)
    output$call = this.call

    return(output)
}
