compute_lambda.default <- function(PeakO,w1,h1,X,w2,h2,D,
                                   alpha = 0.5,
                                   beta = 0.001,
                                   eps = 1e-10,
                                   ...) {
    return(computeLambda(PeakO,w1,h1,X,w2,h2,D,
                         alpha,
                         beta,
                         eps))
}
