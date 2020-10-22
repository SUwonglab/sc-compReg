#include "../inst/include/coupled_nmf.h"


// [[Rcpp::export]]
Rcpp::List initializeMatrix(const unsigned int POnRow,
                                 const unsigned int POnCol,
                                 const unsigned int XnCol,
                                 const unsigned int k,
                                 const arma::sp_mat& D) {
    try {
        arma::mat W10 = arma::mat(POnRow, k, fill::randu);
        arma::mat H10 = arma::mat(k, POnCol, fill::randu);

        arma::mat W20 = D * W10;
        W20 = W20 - W20.min();
        arma::mat H20 = arma::mat(k, XnCol, fill::randu);

        arma::rowvec W10Sum = arma::sum(arma::pow(W10, 2), 0); // dim=0 to sum over columns
        arma::rowvec W20Sum = arma::sum(arma::pow(W20, 2), 0);
        arma::mat H1 = arma::diagmat(W10Sum) * H10;
        arma::mat H2 = arma::diagmat(W20Sum) * H20;
        arma::mat W1 = W10 * arma::diagmat(1./ arma::sqrt(W10Sum));
        arma::mat W2 = W20 * arma::diagmat(1./ arma::sqrt(W20Sum));

        return Rcpp::List::create(Named("W1") = W1,
                                 Named("W2") = W2,
                                 Named("H1") = H1,
                                 Named("H2") = H2);
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}


// [[Rcpp::export]]
Rcpp::List computeLambda(const arma::sp_mat& PeakO,
                         const arma::mat& w1,
                         const arma::mat& h1,
                         const arma::sp_mat& X,
                         const arma::mat& w2,
                         const arma::mat& h2,
                         const arma::sp_mat& D,
                         double alpha,
                         double beta,
                         double eps) {
    try {
        // take mean of both
        double r1 = arma::mean(arma::mean(PeakO * h1.t())) / arma::mean(arma::mean(D.t() * w2));
        double lambda2 = 2 * alpha * beta * r1;

        double r2 = arma::mean(arma::mean(X * h2.t())) / arma::mean(arma::mean(D * w1));
        double lambda1 = alpha / (1 - alpha + eps) * r1 / (r2 + eps);
        return Rcpp::List::create(Named("lambda1") = lambda1,
                                  Named("lambda2") = lambda2);
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}

arma::mat quantileNorm(arma::mat input) {
    // sorting each column
    arma::mat sortedInput = arma::sort(input, "ascend", 0);
    // average over rows
    unsigned int nrows = input.n_rows;
    arma::vec rank = arma::mean(sortedInput, 1);
    arma::mat output = arma::mat(input.n_rows, input.n_cols, arma::fill::zeros);
    arma::uvec index = arma::regspace<uvec>(0, input.n_rows-1);
    arma::vec col;
    unsigned int ind;
    // efficient computing for sparse matrix
    arma::uvec positiveInd = arma::find(input > 0);
    unsigned int c, r;
    for (unsigned int i = 0; i < positiveInd.n_elem; ++i) {
        ind = positiveInd(i);
        if (c != ind / nrows) {
            c = ind / nrows;
            col = sortedInput.col(c);
        }
        r = ind % nrows;
        output(r, c) = rank(arma::min(index.elem(arma::find(col == input(r,c)))));
    }
    return output;
}


// [[Rcpp::export]]
Rcpp::List iterateCluster(const arma::sp_mat& PeakO,
                          const arma::sp_mat& X,
                          const arma::sp_mat& D,
                          const unsigned int k,
                          const unsigned int maxIter,
                          double lambda1,
                          double lambda2,
                          arma::mat W10,
                          arma::mat H10,
                          arma::mat W20,
                          arma::mat H20,
                          double tolX,
                          double tolFun,
                          bool verbose,
                          int loopUpdate
) {
    try {
        double sqrtEps = sqrt(arma::datum::eps);
        unsigned int n = PeakO.n_rows;
        unsigned int m = PeakO.n_cols;
        unsigned int n1 = X.n_rows;
        unsigned int m1 = X.n_cols;

        arma::mat s = arma::eye(k, k);

        double dnorm0 =
                pow(arma::norm(PeakO - W10 * H10, "fro"), 2.0) + lambda1 * pow(arma::norm(X - W20 * H20, "fro"), 2.0);
        double dnorm, dw1, dh1, dw2, dh2, delta;
        arma::mat S1, numer, H1, mu11, mu12, S2, H2, W1, W2, mu21, mu22;
        arma::urowvec S20, S10;
        arma::mat FC1 = arma::mat(PeakO.n_rows, k, arma::fill::zeros);
        arma::mat FC2 = arma::mat(X.n_rows, k, arma::fill::zeros);
        arma::vec tempSumS10 = arma::vec(PeakO.n_rows);
        arma::vec tempSumS20 = arma::vec(X.n_rows);
        arma::vec tempSumNS10 = arma::vec(PeakO.n_rows);
        arma::vec tempSumNS20 = arma::vec(X.n_rows);
        unsigned int SSize, XSize;
        arma::sp_mat zeroPOVec = arma::sp_mat(PeakO.n_rows, 1);
        arma::sp_mat zeroXVec = arma::sp_mat(X.n_rows, 1);
        arma::uvec assignment;
        arma::mat S, WP1, WP2;
        arma::uvec tempJIdx, tempNJIdx;
        arma::vec hTempVec;
        arma::mat H1T, H2T;
        for (int iter = 1; iter <= maxIter; ++iter) {
            if (verbose) std::cout<< "On iteration " << iter << "." << std::endl;
            Rcpp::checkUserInterrupt();
            S1 = 0.5 * lambda2 * D.t() * W20 * s;
            numer = W10.t() * PeakO;
            H1 = H10 % (numer / ((W10.t() * W10) * H10 + arma::eps(numer) + arma::datum::eps));
            H1.elem(arma::find(H1 < 0)).fill(0.);
            numer = PeakO * H1.t() + S1;
            numer.elem(arma::find(numer < 0.)).fill(0.);
            mu11 = arma::diagmat(arma::sum(0.5 * W10 % numer));
            mu12 = arma::diagmat(arma::sum(0.5 * W10 % (W10 * (H1 * H1.t()))));
            W1 = W10 % ((numer + 2 * W10 * mu12) /
                        (W10 * (H1 * H1.t()) + 2 * W10 * mu11 + arma::eps(numer) + arma::datum::eps));

            W1.elem(arma::find(W1 < 0.)).fill(0.);
            S2 = 0.5 * (lambda2 / (lambda1 + arma::datum::eps)) * (D * W1 * s);
            numer = W20.t() * X;
            H2 = H20 % (numer / ((W20.t() * W20) * H20 + arma::eps(numer) + arma::datum::eps));
            H2.elem(arma::find(H2 < 0.)).fill(0.);
            numer = X * H2.t() + S2;
            numer.elem(arma::find(numer<0.)).fill(0.);
            mu21 = arma::diagmat(arma::sum(0.5 * W20 % numer));
            mu22 = arma::diagmat(arma::sum(0.5 * W20 % (W20 * (H2 * H2.t()))));
            W2 = W20 % ((numer + 2 * W20 * mu22) /
                 (W20 * (H2 * H2.t()) + 2 * W20 * mu21 + arma::eps(numer) + arma::datum::eps));
            W2.elem(arma::find(W2 < 0.)).fill(0.);

            dnorm = pow(arma::norm(PeakO - W1 * H1, "fro"), 2.0) + lambda1 * pow(arma::norm(X - W2 * H2, "fro"), 2.0);
            dw1 = arma::max(arma::max(arma::abs(W1 - W10))) / (sqrtEps + arma::max(arma::max(arma::abs(W10))));
            dh1 = arma::max(arma::max(arma::abs(H1 - H10))) / (sqrtEps + arma::max(arma::max(arma::abs(H10))));
            dw2 = arma::max(arma::max(arma::abs(W2 - W20))) / (sqrtEps + arma::max(arma::max(arma::abs(W20))));
            dh2 = arma::max(arma::max(arma::abs(H2 - H20))) / (sqrtEps + arma::max(arma::max(arma::abs(H20))));
            delta = std::max({dw1, dh1, dw2, dh2});

            if (iter > 1) {
                // for debugging purposes for now
                if (delta <= tolX) {
                    std::cout << std::fixed; // print out 6 decimal places
                    std::cout << "delta " << delta << " is small" << endl;
                    break;
                } else if (abs(dnorm0 - dnorm) <= tolFun * std::max(1.0, dnorm0)) {
                    std::cout << std::fixed; // print out 6 decimal places
                    std::cout << "dnorm0 - dnorm " << abs(dnorm0 - dnorm) << " is small" << endl;
                    break;
                } else if (iter == maxIter) {
                    break;
                }
            }

            W10 = W1;
            H10 = H1;
            W20 = W2;
            H20 = H2;

            if (iter % loopUpdate == 0) {
                hTempVec = 1. / arma::sqrt(arma::sum(arma::pow(H20, 2.0), 1));
                H2T = H20;
                // row-wise scalar multiplication
                for (unsigned int hIdx = 0; hIdx < H20.n_rows; ++hIdx) {
                    H2T.row(hIdx) *= hTempVec.at(hIdx);
                }
                S20 = arma::index_max(H2T);
                hTempVec = 1. / arma::sqrt(arma::sum(arma::pow(H10, 2.0), 1));
                H1T = H10;
                for (unsigned int hIdx = 0; hIdx < H10.n_rows; ++hIdx) {
                    H1T.row(hIdx) *= hTempVec.at(hIdx);
                }
                S10 = arma::index_max(H1T);
                SSize = S10.n_elem;
                XSize = S20.n_elem;
                FC1.fill(arma::fill::zeros); // could remove later
                FC2.fill(arma::fill::zeros);
                for (unsigned int j = 0; j < k; ++j) {
                    tempSumS10.fill(arma::fill::zeros);
                    tempSumNS10.fill(arma::fill::zeros);
                    tempSumS20.fill(arma::fill::zeros);
                    tempSumNS20.fill(arma::fill::zeros);

                    for (unsigned int poColIdx = 0; poColIdx < SSize; ++poColIdx) {
                        zeroPOVec = PeakO.col(poColIdx);
                        if (S10.at(poColIdx) == j) {
                            tempSumS10 += zeroPOVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                        } else {
                            tempSumNS10 += zeroPOVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                        }
                    }

                    for (unsigned int xColIdx = 0; xColIdx < XSize; ++xColIdx) {
                        zeroXVec = X.col(xColIdx);
                        if (S20.at(xColIdx) == j) {
                            tempSumS20 += zeroXVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                        } else {
                            tempSumNS20 += zeroXVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                        }
                    }
                    FC1.col(j) = tempSumS10 / (tempSumNS10 / arma::sum(S10 != j) * arma::sum(S10 == j) + 1.0);
                    FC2.col(j) = tempSumS20 / (tempSumNS20 / arma::sum(S20 != j) * arma::sum(S20 == j) + 1.0);
                }
                WP1 = quantileNorm(FC1);
                WP2 = quantileNorm(FC2);
                S = (WP2.t() * D * WP1).t();
                assignment = arma::find(hungarian(-S) == 1);
                assignment = assignment - arma::floor(assignment / S.n_cols) * S.n_cols;
                W2 = W2.cols(assignment);
                H2 = H2.rows(assignment);
            }
            dnorm0 = dnorm;
        }
        W1 = W10;
        H1 = H10;
        W2 = W20;
        H2 = H20;
        double score =
                pow(arma::norm(PeakO - W1 * H1, "fro"), 2.0) + lambda1 * pow(arma::norm(X - W2 * H2, "fro"), 2.0) -
                lambda2 * arma::trace(W2.t() * D * W1);

        H1T = H1;
        H2T = H2;
        hTempVec = 1. / arma::sqrt(arma::sum(arma::pow(H2, 2.0), 1));
        for (unsigned int hIdx = 0; hIdx < H2.n_rows; ++hIdx) {
            H2T.row(hIdx) *= hTempVec.at(hIdx);
        }
        S20 = arma::index_max(H2T);

        hTempVec = 1. / arma::sqrt(arma::sum(arma::pow(H1, 2.0), 1));
        for (unsigned int hIdx = 0; hIdx < H1.n_rows; ++hIdx) {
            H1T.row(hIdx) *= hTempVec.at(hIdx);
        }
        S10 = arma::index_max(H1T);

        FC1.fill(arma::fill::zeros);
        FC2.fill(arma::fill::zeros);
        SSize = S10.n_elem;
        XSize = S20.n_elem;
        for (int j = 0; j < k; ++j) {
            tempSumS10.fill(arma::fill::zeros);
            tempSumNS10.fill(arma::fill::zeros);
            tempSumS20.fill(arma::fill::zeros);
            tempSumNS20.fill(arma::fill::zeros);

            for (unsigned int poColIdx = 0; poColIdx < SSize; ++poColIdx) {
                zeroPOVec = PeakO.col(poColIdx);
                if (S10.at(poColIdx) == j) {
                    tempSumS10 += zeroPOVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                } else {
                    tempSumNS10 += zeroPOVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                }
            }

            for (unsigned int xColIdx = 0; xColIdx < XSize; ++xColIdx) {
                zeroXVec = X.col(xColIdx);
                if (S20.at(xColIdx) == j) {
                    tempSumS20 += zeroXVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                } else {
                    tempSumNS20 += zeroXVec.transform([](double val) { return (val > 0) ? 1. : 0.; });
                }
            }

            FC1.col(j) = tempSumS10 / (tempSumNS10 / arma::sum(S10 != j) * arma::sum(S10 == j) + 1.0);
            FC2.col(j) = tempSumS20 / (tempSumNS20 / arma::sum(S20 != j) * arma::sum(S20 == j) + 1.0);
        }
        WP1 = quantileNorm(FC1);
        WP2 = quantileNorm(FC2);
        S = (WP2.t() * D * WP1).t();
        assignment = arma::find(hungarian(-S) == 1);
        assignment = assignment - arma::floor(assignment / S.n_cols) * S.n_cols;
        W2 = W2.cols(assignment);
        H2 = H2.rows(assignment);
        return Rcpp::List::create(Named("W1") = W1,
                                  Named("W2") = W2,
                                  Named("H1") = H1,
                                  Named("H2") = H2,
                                  Named("score") = score);
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}

// [[Rcpp::export]]
Rcpp::List cluster(const arma::mat& H1,
                   const arma::mat& H2) {
    arma::urowvec C1 = arma::index_max(H1, 0);
    arma::urowvec C2 = arma::index_max(H2, 0);
    return Rcpp::List::create(Named("c1") = C1,
                              Named("c2") = C2);
}


// [[Rcpp::export]]
Rcpp::List postLapMatMult(arma::mat W1,
                          arma::mat W2,
                          arma::mat H1,
                          arma::mat H2
        ) {
    try {
        H1 = arma::diagmat(1. / arma::sqrt(arma::sum(arma::pow(H1, 2.0), 1))) * H1;
        H2 = arma::diagmat(1. / arma::sqrt(arma::sum(arma::pow(H2, 2.0), 1))) * H2;
        W1 = W1 * arma::diagmat(1. / arma::sqrt(arma::sum(arma::pow(W1, 2.0))));
        W2 = W2 * arma::diagmat(1. / arma::sqrt(arma::sum(arma::pow(W2, 2.0))));
        return Rcpp::List::create(Named("W1") = W1,
                                  Named("W2") = W2,
                                  Named("H1") = H1,
                                  Named("H2") = H2);
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}

