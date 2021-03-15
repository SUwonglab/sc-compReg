/** This implementation is from https://gallery.rcpp.org/articles/minimal-assignment/
 * by Lars Simon Zehnder.
 */

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <iostream>
#include <vector>
#include<RcppArmadillo.h>

using namespace std;

void step_one(unsigned int &step, arma::mat &cost,
              const unsigned int &N);

void step_two (unsigned int &step, const arma::mat &cost,
               arma::umat &indM, arma::ivec &rcov,
               arma::ivec &ccov, const unsigned int &N);

void step_three(unsigned int &step, const arma::umat &indM,
                arma::ivec &ccov, const unsigned int &N);

void find_noncovered_zero(int &row, int &col,
                         const arma::mat &cost, const arma::ivec &rcov,
                         const arma::ivec &ccov, const unsigned int &N);

bool star_in_row(int &row, const arma::umat &indM,
                 const unsigned int &N);

void find_star_in_row (const int &row, int &col,
                       const arma::umat &indM, const unsigned int &N);

void step_four (unsigned int &step, const arma::mat &cost,
                arma::umat &indM, arma::ivec &rcov, arma::ivec &ccov,
                int &rpath_0, int &cpath_0, const unsigned int &N);

void find_star_in_col (const int &col, int &row,
                       const arma::umat &indM, const unsigned int &N);

void find_prime_in_row (const int &row, int &col,
                        const arma::umat &indM, const unsigned int &N);

void augment_path (const int &path_count, arma::umat &indM,
                   const arma::imat &path);

void clear_covers (arma::ivec &rcov, arma::ivec &ccov);

void erase_primes(arma::umat &indM, const unsigned int &N);

void step_five (unsigned int &step,
                arma::umat &indM, arma::ivec &rcov,
                arma::ivec &ccov, arma::imat &path,
                int &rpath_0, int &cpath_0,
                const unsigned int &N);

void find_smallest (double &minval, const arma::mat &cost,
                    const arma::ivec &rcov, const arma::ivec &ccov,
                    const unsigned int &N);

void step_six (unsigned int &step, arma::mat &cost,
               const arma::ivec &rcov, const arma::ivec &ccov,
               const unsigned int &N);

arma::umat hungarian(const arma::mat &input_cost);


#endif
