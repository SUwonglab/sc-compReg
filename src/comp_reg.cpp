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

/**
 * Extend division reminder to vectors
 *
 * @param   a       Dividend
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n)
{
    return a - floor(a/n)*n;
}

void parseMotifTarget(std::string filePath,
                      std::vector<std::string>& stringVec1,
                      std::vector<std::string>& stringVec2,
                      std::vector<float>& doubleVec) {
    char buffer[BUFFER_SIZE];
    char a[100], b[100];
    float c;
    int scanRet;
    FILE* f = fopen(filePath.c_str(), "r");
    while (true) {
        if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
        scanRet = sscanf(buffer, "%s %s %lf", a, b, &c);
        stringVec1.push_back(a);
        stringVec2.push_back(b);
        doubleVec.push_back(c);
    }
    fclose(f);
}

// [[Rcpp::export]]
Rcpp::List mfbs(std::vector<std::string> TFName,
                std::vector<std::string> motifName,
                arma::vec motifWeight,
                std::vector<std::string> elementName,
                std::vector<std::string> match2TF,
                std::vector<std::string> match2Motif,
                std::string motifTargetPath) {
    try {
        std::vector<std::string> strVec1, strVec2;
        std::vector<float> f3;
        parseMotifTarget(motifTargetPath, strVec1, strVec2, f3);




    } catch(...) {
        ::Rf_error("c++ exception");
    }
}






// [[Rcpp::export]]
Rcpp::List subpopulationLink(arma::mat EMH,
                             arma::mat EMC,
                             arma::mat OMH,
                             arma::mat OMC) {
    // EMH - EMeanHealthy
    // EMC - EMeanCll
    // OMH - OMeanHealthy
    // OMC - OMeanCll
    try {
        // r1 is K1 x K2
        arma::mat r1 = arma::cor(EMH, EMC);
        arma::mat r2 = arma::cor(EMH, EMC);
        // outer product
        arma::mat rr1 = r1 - arma::sum(r1, 0).t() * arma::sum(r1, 1) / arma::accu(r1);
        arma::mat rr2 = r2 - arma::sum(r2, 0).t() * arma::sum(r2, 1) / arma::accu(r2);
        arma::mat rr = rr1 + rr2;
        arma::uvec b = arma::find(rr > 0);

        arma::mat a;
        arma::uvec rrPos = arma::find(rr > 0);
        unsigned int rrNRows = rr.n_rows;
        // row
        a.col(0) = conv_to<vec>::from(rrPos / rrNRows);
        // col
        a.col(1) = conv_to<vec>::from(mod(rrPos, rrNRows));
        a.col(2) = r1.elem(b);
        a.col(3) = r2.elem(b);
        a.col(4) = rr1.elem(b);
        a.col(5) = rr2.elem(b);
        a.col(6) = rr.elem(b);
        arma::uvec f = arma::sort_index(a.col(6), "descend");
        a = a.rows(f);
        arma::uvec aRows =  arma::find(((a.col(4) > 0) % (a.col(5) > 0)) == 1);
        a = a.rows(aRows);

        unsigned int aNRows = a.n_rows;
        arma::vec S1 = arma::vec(aNRows, arma::fill::zeros);
        arma::vec S2 = arma::vec(aNRows, arma::fill::zeros);
        arma::mat match = arma::mat(aNRows, a.n_cols, arma::fill::zeros);
        unsigned int idx;
        unsigned int matchIdx;
        for (int i = 0; i < aNRows; ++ i) {
            idx = arma::any(S1 == a.at(i, 1)) + arma::any(S2 == a.at(i, 2));
            S1.at(i) = a.at(i, 1);
            S2.at(i) = a.at(i, 2);
            if (idx < 2) {
                match.row(matchIdx) = a.row(i);
                ++matchIdx;
            }
        }
        match.shed_rows(matchIdx, aNRows - 1);
        return Rcpp::List::create(Named("match") = match);
    } catch (...) {
        ::Rf_error("c++ exception");
    }
}



// [[Rcpp::export]]
Rcpp::List clusterProfile(const arma::sp_mat& O1,
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
                   std::vector<std::vector<std::string>> peakNameIntersect) {
    try {
        std::vector<std::string> peakNameIntersect1 = peakNameIntersect.at(0);
        std::vector<std::string> peakNameIntersect2 = peakNameIntersect.at(1);

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
        unsigned int OMeanNRows = m + d1Zero + d2Zero;

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

        for (arma::uvec::iterator it = d1ZeroIdx.begin(); it != d1ZeroIdx.end(); ++it) {
            peakNameIntersect1.push_back(peakName1.at(*it));
        }

        for (arma::uvec::iterator it = d2ZeroIdx.begin(); it != d2ZeroIdx.end(); ++it) {
            peakNameIntersect1.push_back(peakName2.at(*it));
        }

        return Rcpp::List::create(Named("E1Mean") = E1Mean,
                                  Named("E2Mean") = E2Mean,
                                  Named("symbol") = symbol,
                                  Named("O1Mean") = O1Mean,
                                  Named("O2Mean") = O2Mean,
                                  Named("elementName") = peakNameIntersect1);
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

