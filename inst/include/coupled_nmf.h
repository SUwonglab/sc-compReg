//coupledNMF.h
#ifndef COUPLED_NMF
#define COUPLED_NMF

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <random>
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath> // isnan
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <math.h>
#include <chrono> // for seeding
#include "Hungarian.h"

#include <iomanip>

using namespace RcppParallel;
using namespace arma;
using namespace R;
using namespace Rcpp;

#endif //COUPLED_NMF
