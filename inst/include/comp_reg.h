//comp_reg.h
#ifndef COMP_REG
#define COMP_REG

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath> // isnan
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <math.h>
#include <chrono> // for seeding
#include <bits/stdc++.h>
#include <cstdio>
#include <algorithm>
#include <tuple>
#include <boost/math/distributions/students_t.hpp> // t-test
#include <boost/math/distributions/chi_squared.hpp> // LRT
#include <boost/math/distributions/gamma.hpp> // gamma quantile matching

#include <iomanip>

#define BUFFER_SIZE 100
#define TOLERANCE 0.00000001
#define ALPHA_THRESH 0.05

using namespace RcppParallel;
using namespace arma;
using namespace R;
using namespace Rcpp;

#endif //COMP_REG
