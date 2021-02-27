mfbs.default <- function(elem.name,
                         symbol,
                         motif.name,
                         motif.weight,
                         match2,
                         motif.file) {
    this.call <- match.call()
    tf.name <- intersect(symbol, unique(match2[[2]]))
    tf.name.len <- length(tf.name)
    tf.binding <- matrix(0, tf.name.len, length(elem.name))

    f3 <- motif.file$C3

    f1 <- match(motif.file$C1, elem.name, nomatch = 0)
    d1 <- f1 > 0

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


    motif.binding <- sparseMatrix(dims = c(length(motif.name),length(elem.name)),
                                  i = f2,
                                  j = f1,
                                  x = f3)
    motif.weight.len <- length(motif.weight)


    motif.binding <- mult(Matrix(diag(1 / (motif.weight + 0.1),
                                      nrow = motif.weight.len,
                                      ncol = motif.weight.len), sparse=T), motif.binding)
    motif.binding@x <- log(motif.binding@x)

    mf1 <- match(match2[[1]], motif.name, nomatch=0)
    mf2 <- match(match2[[2]], tf.name, nomatch=0)

    for (i in 1:tf.name.len) {
        a <- which(mf2 == i)
        if (length(a) > 1) {
            tf.binding[i, ] <- colMax(motif.binding[mf1[a], ])
        } else if (length(a) == 1) {
            tf.binding[i, ] <- motif.binding[mf1[a], ]
        }
    }

    output <- list()
    output$tf.binding <- tf.binding
    output$tf.name <- tf.name
    output$call = this.call
    return(output)
}
