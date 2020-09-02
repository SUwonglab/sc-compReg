#include "../inst/include/comp_reg.h"


std::vector<std::string> intersection(std::vector<std::string> &v1,
                                      std::vector<std::string> &v2){
    std::vector<std::string> v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}


// [[Rcpp::export]]
Rcpp::List compReg(const arma::sp_mat& O1,
                   const arma::sp_mat& E1,
                   arma::uvec O1Idx,
                   arma::uvec E1Idx,
                   std::vector<std::string> symbol1,
                   std::vector<std::string> peakName1,
                   const arma::sp_mat& O2,
                   const arma::sp_mat& E2,
                   arma::uvec O2Idx,
                   arma::uvec E2Idx,
                   std::vector<std::string> symbol2,
                   std::vector<std::string> peakName2,
                   std::vector<std::string> peakNameIntersect1,
                   std::vector<std::string> peakNameIntersect2
                   ) {
    try {
        unsigned int K1 = std::max(arma::max(O1Idx), arma::max(E1Idx));
        unsigned int K2 = std::max(arma::max(O2Idx), arma::max(E2Idx));
        auto symbol = intersection(symbol1, symbol2);
        std::sort(symbol1.begin(), symbol1.end());
        std::vector<std::string>::iterator s1Start= symbol1.begin();
        std::sort(symbol2.begin(), symbol2.end());
        std::vector<std::string>::iterator s2Start= symbol2.begin();

        arma::uvec f2 = arma::uvec(symbol.size(), arma::fill::zeros);
        arma::uvec f1 = arma::uvec(symbol.size(), arma::fill::zeros);
        unsigned int ind = 0;
        for (std::vector<std::string>::iterator t=symbol.begin(); t!=symbol.end(); ++t) {
            f1(ind) = lower_bound(symbol1.begin(), symbol1.end(), *t) - s1Start;
            f2(ind) = lower_bound(symbol2.begin(), symbol2.end(), *t) - s2Start;
            ++ind;
        }

        arma::mat E1Mean = arma::mat(f1.n_rows, K1, arma::fill::zeros);
        arma::mat E2Mean = arma::mat(f2.n_rows, K2, arma::fill::zeros);

        double avg = 0.;
        unsigned int f1Idx;
        for (int i = 0; i < K1; ++i) {
            arma::uvec idx = arma::find(E1Idx == i);
            f1Idx = 0;
            for (arma::uvec::iterator it = f1.begin(); it != f1.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    avg += E1.row(*it).at(*elemIt);
                }
                E1Mean.at(f1Idx, i) = avg / idx.n_elem;
                ++f1Idx;
            }
        }

        unsigned int f2Idx;
        for (int i = 0; i < K2; ++i) {
            arma::uvec idx = arma::find(E2Idx == i);
            f2Idx = 0;
            for (arma::uvec::iterator it = f2.begin(); it != f2.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    avg += E2.row(*it).at(*elemIt);
                }
                E2Mean.at(f1Idx, i) = avg / idx.n_elem;
                ++f2Idx;
            }
        }

        std::sort(peakName1.begin(), peakName1.end());
        std::vector<std::string>::iterator p1Start= peakName1.begin();
        std::sort(peakName2.begin(), peakName2.end());
        std::vector<std::string>::iterator p2Start= peakName2.begin();

        f2 = arma::uvec(peakNameIntersect2.size(), arma::fill::zeros);
        f1 = arma::uvec(peakNameIntersect1.size(), arma::fill::zeros);
        ind = 0;
        for (std::vector<std::string>::iterator t=peakNameIntersect1.begin(); t!=peakNameIntersect1.end(); ++t) {
            f1(ind) = lower_bound(peakName1.begin(), peakName1.end(), *t) - s1Start;
            ++ind;
        }
        ind = 0;
        for (std::vector<std::string>::iterator t=peakNameIntersect2.begin(); t!=peakNameIntersect2.end(); ++t) {
            f2(ind) = lower_bound(peakName2.begin(), peakName2.end(), *t) - s2Start;
            ++ind;
        }
        std::sort(peakNameIntersect1.begin(), peakNameIntersect1.end());
        std::sort(peakNameIntersect2.begin(), peakNameIntersect2.end());
        arma::uvec d1 = arma::uvec(peakName1.size(), arma::fill::zeros);
        arma::uvec d2 = arma::uvec(peakName2.size(), arma::fill::zeros);
        ind = 0;
        for (std::vector<std::string>::iterator t=peakName1.begin(); t!=peakName1.end(); ++t) {
            d1.at(ind) = std::binary_search(peakNameIntersect1.begin(), peakNameIntersect1.end(), *t);
            ++ind;
        }
        ind = 0;
        for (std::vector<std::string>::iterator t=peakName2.begin(); t!=peakName2.end(); ++t) {
            d2.at(ind) = std::binary_search(peakNameIntersect2.begin(), peakNameIntersect2.end(), *t);
            ++ind;
        }

        unsigned int d1Zero, d2Zero;
        arma::uvec d1ZeroIdx = arma::find(d1 == 0);
        d1Zero = d1ZeroIdx.n_elem;
        arma::uvec d2ZeroIdx = arma::find(d2 == 0);
        d2Zero = d2ZeroIdx.n_elem;
        unsigned int m = peakNameIntersect1.size(); // length(peakNameIntersect1) = length(peakNameIntersect2)
        const std::vector<unsigned int> v {peakNameIntersect1.size(), d1Zero, d2Zero};
        auto maxElem = std::max_element(v.begin(), v.end());
        unsigned int OMeanNRows = *maxElem; // TODO: check with Duren here

        arma::mat O1Mean = arma::mat(OMeanNRows, K1, arma::fill::zeros);
        avg = 0.;
        f1Idx = 0;
        double loopElem;
        for (int i = 0; i < K1; ++i) {
            arma::uvec idx = arma::find(O1Idx == i);
            f1Idx = 0;
            for (arma::uvec::iterator it = f1.begin(); it != f1.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    loopElem = O1.row(*it).at(*elemIt);
                    if (loopElem > 0) avg += loopElem;
                }
                if (f1Idx < m) O1Mean.at(f1Idx, i) = avg / idx.n_elem;
                ++f1Idx;
            }

            f1Idx = m;
            for (arma::uvec::iterator it = d1ZeroIdx.begin(); it != d1ZeroIdx.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    loopElem = O1.row(*it).at(*elemIt);
                    if (loopElem > 0) avg += loopElem;
                }
                O1Mean.at(f1Idx, i) = avg / idx.n_elem;
                ++f1Idx;
            }
        }

        arma::mat O2Mean = arma::mat(OMeanNRows, K2, arma::fill::zeros);
        avg = 0.;
        f2Idx = 0;
        loopElem;
        for (int i = 0; i < K2; ++i) {
            arma::uvec idx = arma::find(O2Idx == i);
            f2Idx = 0;
            for (arma::uvec::iterator it = f2.begin(); it != f2.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    loopElem = O2.row(*it).at(*elemIt);
                    if (loopElem > 0) avg += loopElem;
                }
                if (f2Idx < m) O2Mean.at(f2Idx, i) = avg / idx.n_elem;
                ++f2Idx;
            }

            f2Idx = m;
            for (arma::uvec::iterator it = d1ZeroIdx.begin(); it != d1ZeroIdx.end(); ++it) {
                avg = 0;
                for (arma::uvec::iterator elemIt = idx.begin(); elemIt != idx.end(); ++elemIt) {
                    loopElem = O2.row(*it).at(*elemIt);
                    if (loopElem > 0) avg += loopElem;
                }
                O2Mean.at(f2Idx, i) = avg / idx.n_elem;
                ++f1Idx;
            }
        }

        O1Mean /= arma::mean(O1Mean, 0);
        O2Mean /= arma::mean(O2Mean, 0);

        return Rcpp::List::create(Named("E1Mean") = E1Mean,
                                  Named("E2Mean") = E2Mean,
                                  Named("symbol") = symbol,
                                  Named("O1Mean") = O1Mean,
                                  Named("O2Mean") = O2Mean,
                                  Named("d1Zero") = d1ZeroIdx + 1,
                                  Named("d2Zero") = d2ZeroIdx + 1);
    } catch(...) {
            ::Rf_error("c++ exception");
    }
}



std::vector<std::string> parseSymbol(std::string filePath) {
    std::ifstream symbolFile(filePath);
    std::vector <std::string> symbols;

    if ( !symbolFile ) {
        std::cerr << "File does not exist!\n";
    }

    std::string line;
    while (getline(symbolFile, line)) {
        if (std::find(symbols.begin(), symbols.end(), line) == symbols.end()) {
            // new symbol is not in file
            symbols.push_back(line);
        }
    }
    symbolFile.close();
    return symbols;
}
