//quantile_norm.h
#ifndef QUANTILE_NORM
#define QUANTILE_NORM

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace R;
using namespace Rcpp;


arma::mat quantileNorm(arma::mat input);

#endif //QUANTILE_NORM
