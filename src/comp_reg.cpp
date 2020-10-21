#include "../inst/include/comp_reg.h"

// [[Rcpp::export]]
Rcpp::List loadPeakNameIntersectFile(std::string path,
                                     char token) {
    try {
        char buffer[BUFFER_SIZE];
        char *a = (char *)malloc(BUFFER_SIZE);
        char *b = (char *)malloc(BUFFER_SIZE);
        std::string tokenStr = std::string(1, token);
        tokenStr = "%s" + tokenStr + "%s";
        std::vector<std::string> file1;
        std::vector<std::string> file2;
        int scanRet;
        FILE* f = fopen(path.c_str(), "r");
        std::string temp;
        while (true) {
            if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
            scanRet = sscanf(buffer, tokenStr.c_str(), a, b);
            temp = a;
            file1.push_back(temp);
            temp = b;
            file2.push_back(temp);
        }
        fclose(f);
        free(a);
        free(b);
        return Rcpp::List::create(Named("vo") = file1,
                                  Named("vt") = file2);
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}

std::vector<std::string> intersection(std::vector<std::string> v1,
                                      std::vector<std::string> v2){
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
    char *a = (char *)malloc(BUFFER_SIZE);
    char *b = (char *)malloc(BUFFER_SIZE);
    float c;
    int scanRet;
    FILE* f = fopen(filePath.c_str(), "r");
    std::string temp;
    while (true) {
        if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
        scanRet = sscanf(buffer, "%s %s %f", a, b, &c);
        temp = a;
        stringVec1.push_back(temp);
        temp = b;
        stringVec2.push_back(temp);
        doubleVec.push_back(c);
    }
    fclose(f);
    free(a);
    free(b);
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
        unsigned int elemNameSize = elementName.size();
        std::vector<std::string> strVec1, strVec2;
        std::vector<float> floatVec3;
        unsigned int index;
        parseMotifTarget(motifTargetPath, strVec1, strVec2, floatVec3);

        // check if and where is every element of strVec1 in elementName
        auto elemEnd = elementName.end();
        auto elemStart = elementName.begin();
        arma::vec d1 = arma::vec(strVec1.size(), arma::fill::zeros);
        arma::vec f1 = arma::vec(strVec1.size(), arma::fill::zeros);
        // TODO: first sort then write ismember will be BUGGY!!!
//        std::sort(elemStart, elemEnd);
        unsigned int armaVecIdx = 0;
        arma::vec t2 = arma::regspace(0, elemNameSize - 1);
        for (auto & str : strVec1) {
            auto it = std::find(elemStart, elemEnd, str);
            if (it != elemEnd) {
                index = distance(elemStart, it);
                d1.at(armaVecIdx) = 1;
                f1.at(armaVecIdx) = index;
                if (arma::any(t2 == index)) {
                    t2.shed_rows(arma::find(t2 == index));
                }
            } else {
                f1.at(armaVecIdx) = 0;
                d1.at(armaVecIdx) = 0;
            }
            ++armaVecIdx;
        }

        armaVecIdx = 0;
        auto motifEnd = motifName.end();
        auto motifStart = motifName.begin();
        arma::vec d2 = arma::vec(strVec2.size(), arma::fill::zeros);
        arma::vec f2 = arma::vec(strVec2.size(), arma::fill::zeros);
//        std::sort(motifStart, motifEnd);
        arma::vec t1 = arma::regspace(0, motifName.size()-1); // regspace is inclusive
        for (auto & str : strVec2) {
            auto it = std::find(motifStart, motifEnd, str);
            if (it != motifEnd) {
                index = distance(motifStart, it);
                d2.at(armaVecIdx) = 1;
                f2.at(armaVecIdx) = index;
                if (arma::any(t1 == index)) {
                    t1.shed_rows(arma::find(t1 == index));
                }
            } else {
                d2.at(armaVecIdx) = 0;
                f2.at(armaVecIdx) = 0;
            }
            ++armaVecIdx;
        }

        arma::uvec d1d2BothOnes = arma::find((d1 % d2) == 1);
        f2 = f2.elem(d1d2BothOnes);
        f2 = arma::join_horiz(f2, t1);
        f1 = f1.elem(d1d2BothOnes);
        f1 = arma::join_horiz(f1, arma::vec(t1.n_elem, arma::fill::ones));
        arma::vec f3 = arma::conv_to<arma::vec>::from(floatVec3);
        f3 = arma::join_horiz(f3, arma::vec(t1.n_elem, arma::fill::zeros));

        f1 = arma::join_horiz(f1, t2);
        f2 = arma::join_horiz(f2, arma::vec(t2.n_elem, arma::fill::ones));
        f3 = arma::join_horiz(f3, arma::vec(t2.n_elem, arma::fill::zeros));
        arma::umat spMatLocation = arma::umat(2, f2.n_elem);
        spMatLocation.col(0) = arma::conv_to<arma::uvec>::from(f2);
        spMatLocation.col(1) = arma::conv_to<arma::uvec>::from(f1);
        arma::sp_mat motifBinding = arma::sp_mat(spMatLocation, f3, motifName.size(), elemNameSize);
        motifBinding = arma::diagmat(1 / (motifWeight + 0.1)) * motifBinding;
        motifBinding.transform([](double val) {return log(val + 1);});

        arma::mat TFBinding = arma::mat(TFName.size(), elemNameSize, arma::fill::zeros);

        arma::vec motifIdxVec = arma::vec(match2Motif.size(), arma::fill::zeros);
        armaVecIdx = 0;
        for (auto & str : match2Motif) {
            auto it = std::find(motifStart, motifEnd, str);
            if (it != motifEnd) {
                index = distance(motifStart, it);
                motifIdxVec.at(armaVecIdx) = index;
            }
            ++armaVecIdx;
        }

        auto TFStart = TFName.begin();
        auto TFEnd = TFName.end();
        arma::vec TFIdxVec = arma::vec(match2TF.size(), arma::fill::zeros);
//        std::sort(TFStart, TFEnd);
        armaVecIdx = 0;
        for (auto & str : match2TF) {
            auto it = std::find(TFStart, TFEnd, str);
            if (it != TFEnd) {
                index = distance(TFStart, it);
                TFIdxVec.at(armaVecIdx) = index;
            }
            ++armaVecIdx;
        }
        arma::uvec a;
        arma::mat tempMotifMat;
        unsigned int numAElem;
        unsigned int tempMotifMatIdx;
        for (unsigned int i = 0; i < TFName.size(); ++i) {
            a = arma::find(TFIdxVec == i);
            numAElem = a.n_elem;
            if (numAElem > 1) {
                tempMotifMat = arma::mat(numAElem, elemNameSize);
                tempMotifMatIdx = 0;
                for (uvec::iterator it = a.begin(); it != a.end(); ++it) {
                    tempMotifMat.row(tempMotifMatIdx) = arma::rowvec(motifBinding.row(motifIdxVec.at(*it)));
                }
                TFBinding.row(i) = arma::max(tempMotifMat);
            } else if (numAElem == 1) {
                TFBinding.row(i) = arma::rowvec(motifBinding.row(a.at(0)));
            } else {
                TFBinding.row(i) = arma::rowvec(elemNameSize, arma::fill::zeros);
            }
        }
        return Rcpp::List::create(Named("TF_binding") = arma::sp_mat(TFBinding));
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}


void parsePeakGenePrior(std::string filePath,
                        std::vector<std::string>& stringVec1,
                        std::vector<std::string>& stringVec2,
                        std::vector<float>& doubleVec1,
                        std::vector<float>& doubleVec2) {
    char buffer[BUFFER_SIZE];
    char *a = (char *)malloc(BUFFER_SIZE);
    char *b = (char *)malloc(BUFFER_SIZE);
    float c;
    float d;
    int scanRet;
    FILE* f = fopen(filePath.c_str(), "r");
    std::string temp;
    while (true) {
        if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
        scanRet = sscanf(buffer, "%s %s %f %f", a, b, &c, &d);
        temp = a;
        stringVec1.push_back(temp);
        temp = b;
        stringVec2.push_back(temp);
        doubleVec1.push_back(c);
        doubleVec2.push_back(d);
    }
    fclose(f);
    free(a);
    free(b);
}

template <typename T>
/**
 * ismember(A,B) returns an array containing
 * logical 1 (true) where the data in A is found in B.
 * Elsewhere, the array contains logical 0 (false).
 * The second return value contains the index of the first
 * instance found in array B.
 * @tparam T
 * @tparam T
 * @param vecA
 * @param vecB
 * @return
 */
std::tuple<arma::vec, arma::vec> isMember(const T& vecA,
                                          const T& vecB) {
    auto ABegin = vecA.begin();
    auto AEnd = vecA.end();
    arma::vec lia = arma::vec(vecA.size(), arma::fill::zeros);
    arma::vec locb = arma::vec(vecB.size(), arma::fill::zeros);

    unsigned int insertIdx = 0;
    for (auto t : vecA) {
        auto p = std::find(ABegin, AEnd, t);
        if (p != AEnd) {
            lia.at(insertIdx) = 1.;
            locb.at(insertIdx) = std::distance(p, AEnd);
        } else {
            lia.at(insertIdx) = 0.;
            locb.at(insertIdx) = arma::datum::nan;
        }
        ++insertIdx;
    }
    return std::tuple<arma::vec, arma::vec>{lia, locb};
}

arma::mat uniqueRows(const arma::mat& m) {
    arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);
    for (arma::uword i = 0; i < m.n_rows; i++) {
        for (arma::uword j = i + 1; j < m.n_rows; j++) {
            if (arma::approx_equal(m.row(i), m.row(j), "absdiff", TOLERANCE)) { ulmt(j) = 1; break; }
        }
    }
    return m.rows(find(ulmt == 0));
}

arma::uvec findUniqueRowIdx(const arma::mat& origMat,
                            const arma::mat& uniqueRows) {
    arma::uvec ic = arma::uvec(origMat.n_rows, arma::fill::zeros);
    for (unsigned int oIdx = 0; oIdx < origMat.n_rows; ++oIdx) {
        for (unsigned int uIdx = 0; uIdx < uniqueRows.n_rows; ++uIdx) {
            if (arma::approx_equal(uniqueRows.row(uIdx), origMat.row(oIdx), "absdiff", TOLERANCE)) {
                ic.at(oIdx) = uIdx;
                break;
            }
        }
    }
    return ic;
}


arma::mat accumArrayMin(const arma::uvec& subs,
                        const arma::mat& val) {
    int maxSubs = arma::max(subs);
    arma::vec retVec = arma::vec(maxSubs, arma::fill::zeros);
    arma::uvec idx;
    for (int i = 0; i < maxSubs; ++i) {
        idx = arma::find(subs == i);
        if (idx.n_elem > 0) {
            retVec.at(i) = arma::min(val.elem(idx));
        }
    }
    return retVec;
}

double ttest(const arma::vec& x,
             const arma::vec& y,
             const unsigned int& n,
             const unsigned int& m,
             const unsigned int& df) {
    double numerator = arma::mean(x) - arma::mean(y);
    double denominator = (n-1) * arma::var(x, 0) + (m-1) * arma::var(y, 0);
    denominator = sqrt(denominator / df);
    return numerator / denominator;
}

arma::vec ttest2(const arma::mat& x,
                 const arma::mat& y) {
    assert(x.n_cols == y.n_cols);
    unsigned int numCols = x.n_cols;
    arma::vec retVec = arma::vec(numCols, arma::fill::zeros);
    unsigned int n = x.n_elem;
    unsigned int m = y.n_elem;
    unsigned int df = n + m - 2;
    boost::math::students_t dist(df);
    for (unsigned int i = 0; i < numCols; ++i) {
        retVec.at(i) = boost::math::cdf (
            boost::math::complement(dist, fabs(ttest(x.col(i), y.col(i), n, m, df))));
    }
    return retVec;
}

arma::vec fdrBH(const arma::vec& pvals) {
    arma::vec pSorted = arma::sort(pvals);
    arma::uvec unsortIdx = arma::sort_index(arma::sort_index(pvals));
    unsigned int m = pSorted.n_elem;
    arma::vec mSeq = arma::regspace(0, 1, m-1);
    arma::vec thresh = mSeq * ALPHA_THRESH / m;
    arma::vec wtdP = m * pSorted / mSeq;

    arma::vec adjPVec = arma::vec(m, arma::fill::zeros);
    arma::vec wtdPSorted = arma::sort(wtdP);
    arma::uvec wtdPIdx = arma::sort_index(wtdP);
    unsigned int nextFill = 0;
    for (unsigned int k = 0; k < m; ++k) {
        if (wtdPIdx.at(k) >= nextFill) {
            adjPVec.rows(nextFill, wtdPIdx.at(k)).fill(wtdPSorted.at(k));
            nextFill = wtdPIdx.at(k) + 1;
            if (nextFill >= m) break;
        }
    }
    adjPVec = adjPVec.rows(unsortIdx);
    return adjPVec;
}

arma::sp_mat selectCols(const arma::sp_mat& input,
                        const arma::uvec& idx) {
    arma::sp_mat retMat = arma::sp_mat(input.n_rows, idx.n_elem);
    for (unsigned int i = 0; i < idx.n_elem; ++i) {
        retMat.col(i) = input.col(idx.at(i));
    }
    return retMat;
}

arma::mat corr(const arma::mat& X,
               const arma::mat& Y) {
    arma::mat pMat = arma::mat(X.n_cols, Y.n_cols, arma::fill::zeros);
    double r;
    unsigned int df = X.n_rows + Y.n_rows - 2;
    boost::math::students_t dist(df);
    double tStat;
    for (unsigned int i = 0; i < X.n_cols; ++i) {
        for (unsigned int j = 0; j < Y.n_cols; ++j) {
            r = arma::as_scalar(arma::cor(X.col(i), Y.col(j)));
            tStat = r * sqrt(df) / sqrt(1 - pow(r, 2));
            pMat.at(i, j) = boost::math::cdf (
                    boost::math::complement(dist, fabs(tStat)));
        }
    }
    return pMat;
}

void convertIdxToRowCol(const arma::uvec& idxVec,
                   arma::umat& idxMat,
                   unsigned int nRows) {
    // matrix is interpretted using column-by-column ordering in armadillo
    idxMat.col(0) = mod(idxVec, nRows);
    // col
    idxMat.col(1) = idxVec / nRows;
}

std::tuple<double, double> LRT(double uLogL,
                               double rLogL,
                              unsigned int dof) {
    double testStat = 2 * (uLogL - rLogL);
    boost::math::chi_squared dist(dof);
    double p = boost::math::cdf(dist, fabs(testStat));
    return std::tuple<double, double> {p, testStat};
}


std::tuple<double, double, double, double> bivariateNormalConditionalLR(const arma::mat& X1,
                                                        const arma::mat& X2) {
    double sigmaMinPseudoCount = 0;
    arma::mat XX = join_vert(X1, X2);
    arma::mat beta = arma::pinv(
            arma::join_horiz(
                    arma::vec(XX.n_rows, arma::fill::ones), XX.col(1)) * XX.col(0));
    arma::vec estX = beta.at(0) + beta.at(1) * XX.col(1);
    double condVar = arma::var(XX.col(0) - estX);
    arma::vec S = arma::vec(estX.n_elem);
    S.fill(sqrt(condVar) + sigmaMinPseudoCount);
    double LM0 = arma::as_scalar(arma::sum(arma::log_normpdf(XX.col(0), estX, S)));

    arma::mat beta1 = arma::pinv(arma::join_horiz(arma::vec(X1.n_rows, arma::fill::ones),
                                                  X1.col(1)) * X1.col(0));
    arma::vec estX1 = beta1.at(0) + beta.at(1) * X2.col(1);
    double condVar1 = arma::var(X1.col(0) - estX1);

    arma::mat beta2 = arma::pinv(arma::join_horiz(
            arma::vec(X2.n_rows, arma::fill::ones), X2.col(1) * X2.col(0)
            ));
    arma::vec estX2 = beta2.at(0) + beta2.at(1) * X2.col(1);
    double condVar2 = arma::var(X2.col(0) - estX2);
    S.fill(sqrt(condVar1) + sigmaMinPseudoCount);
    double LM1One = arma::as_scalar(
            arma::sum(arma::log_normpdf(
                    X1.col(0), estX1, S
                    ))
    );
    S.fill(sqrt(condVar2) + sigmaMinPseudoCount);
    double LM1Two = arma::as_scalar(
            arma::sum(arma::log_normpdf(
                    X2.col(0), estX2, S
                    ))
            );
    double LM1 = LM1One + LM1Two;
    std::tuple<double, double> retTup = LRT(LM1, LM0, 3);
    return std::tuple<double, double, double, double>
              {LM0, LM1, std::get<0>(retTup), std::get<1>(retTup)};
}

arma::vec gammaQuantileMatch(arma::vec& X, arma::vec eita) {
    eita = eita.t();
    arma::vec cut = quantile(X, eita * 100);
    return cut;

}


/**
 * Compute indices based on row and column indices.
 */
inline arma::uvec arr2ind(arma::uvec c, arma::uvec r, unsigned int nrow)
{
    return c * nrow + r;
}

/**
 * Assuming indices are sorted, extract the corresponding entries in input
 * array which should be of std::vector< > type.
 * @tparam T
 * @param input
 * @param indices
 * @return
 */
template <class T>
T extractElems(T& input,
               arma::uvec indices) {
    T retVec;
    for (unsigned int i = 0; i < indices.n_elem; ++i) {
        retVec.push_back(input[indices.at(i)]);
    }
    return retVec;
}


// [[Rcpp::export]]
Rcpp::List compReg(arma::mat TFBinding,
                   arma::mat match,
                   const arma::sp_mat& E1,
                   arma::uvec E1Idx,
                   const arma::sp_mat& E2,
                   arma::uvec E2Idx,
                   const arma::mat& O1Mean,
                   const arma::mat& O2Mean,
                   std::vector<std::string> symbol,
                   std::vector<std::string> TFName,
                   std::vector<std::string> elementName,
                   std::string peakGenePriorPath
                   ) {
    try {
        arma::umat TFIdx = arma::sort_index(TFBinding, 2);
        TFIdx = TFIdx.cols(0, 4999); // indices are inclusive
        arma::mat a = arma::mat(TFIdx.n_rows, TFIdx.n_cols, arma::fill::zeros);
        arma::urowvec TFRow;
        arma::rowvec TFBindingRow;
        for (unsigned int rowIdx=0; rowIdx < TFIdx.n_rows; ++rowIdx) {
            TFRow = TFIdx.row(rowIdx);
            TFBindingRow = TFBinding.row(rowIdx);
            a.row(rowIdx) = TFBindingRow.elem(TFRow);
        }
        TFBinding.elem(arma::find(TFBinding - a.col(a.n_cols - 1) < 0)).fill(0);

        std::vector<std::string> strVec1, strVec2;
        std::vector<float> floatVec3, floatVec4;
        parsePeakGenePrior(peakGenePriorPath, strVec1, strVec2, floatVec3, floatVec4);
        arma::vec fv3 = arma::conv_to<arma::vec>::from(floatVec3);
        arma::vec fv4 = arma::conv_to<arma::vec>::from(floatVec4);

//        auto elemBegin = elementName.begin();
//        auto elemEnd = elementName.end();
//
//        arma::vec d = arma::vec(strVec1.size(), arma::fill::zeros);
//        arma::vec f = arma::vec(strVec1.size(), arma::fill::zeros);
//        unsigned int idx = 0;
//        for (auto t : strVec1) {
//            auto p = std::find(elemBegin, elemEnd, t);
//            if (p != elemEnd) {
//                d.at(idx) = 1.;
//                f.at(idx) = std::distance(p, elemEnd);
//            } else {
//                d.at(idx) = 0.;
//                f.at(idx) = NAN;
//            }
//            ++idx;
//        }

        std::tuple<arma::vec, arma::vec> tup = isMember(strVec1, elementName);
        arma::vec d = std::get<0>(tup);
        arma::vec f = std::get<1>(tup);

        std::tuple<arma::vec, arma::vec> tup2 = isMember(strVec2, symbol);
        arma::vec d1 = std::get<0>(tup2);
        arma::vec f1 = std::get<1>(tup2);

        // TODO check whether the .* is outer product
        arma::mat ff = arma::join_horiz(f.elem(arma::find(d % d1 == 1)), f1.elem(arma::find(d % d1 == 1)));
        arma::mat f2 = uniqueRows(ff);
        arma::uvec ic = findUniqueRowIdx(ff, f2);

        arma::vec c3 = accumArrayMin(ic, fv3.elem(arma::find(d % d1 == 1)));
        arma::vec c4 = accumArrayMin(ic, fv4.elem(arma::find(d % d1 == 1)));

        c4.elem(arma::find(c4 < 0.2)).fill(0.);
        arma::vec c = arma::exp(-1. * c3 / 500000.) % c4;
        arma::umat temp = arma::conv_to<arma::umat>::from(f2);
        arma::umat loc = arma::join_cols(temp.col(1), temp.col(0));
        arma::sp_mat beta = arma::sp_mat(loc, c, symbol.size(), elementName.size());

        double i1, i2;
        arma::mat BO1, BO2, TG1, TG2, TF, corrP1, corrP2, pCombine;
        arma::mat OTF1, OTF2, X1, X2;
        unsigned int n1, n2;
        arma::vec pTTest, adjPTTest, LRi, phat, pGamma, LRSColTwo, adjPGamma;
        arma::uvec diffGene, id1, id, tempIdx, sortedIdx;
        arma::umat netIdx;
        std::tuple<double, double, double, double> bivarRetTup;
        arma::mat LRSummary = arma::mat(1, 4, arma::fill::zeros);
        double tenSixteenthPow = (double) pow(10, -16);
        for (unsigned int ii = 0; ii < match.n_rows; ++ii) {
            i1 = match.at(ii, 0);
            i2 = match.at(ii, 2);
            BO1 = (TFBinding % O1Mean.col(i1).t()) * beta.t();
            BO2 = (TFBinding % O2Mean.col(i2).t()) * beta.t();
            TG1 = selectCols(E1, arma::find(E1Idx == i1));
            TG2 = selectCols(E2, arma::find(E2Idx == i2));
            // mean of TG1 and mean of TG2
            TG2 = TG2 * (arma::accu(TG1) / TG1.n_elem) / (arma::accu(TG2) / TG2.n_elem);
            tup = isMember(TFName, symbol);
            d = std::get<0>(tup);
            f = std::get<1>(tup);
            TF = arma::join_horiz(TG1.rows(arma::conv_to<arma::uvec>::from(f)), TG2.rows(arma::conv_to<arma::uvec>::from(f)));
            n1 = TG1.n_cols;
            n2 = TG2.n_cols;
            pTTest = ttest2(TG1.t(), TG2.t());
            adjPTTest = fdrBH(pTTest);
            diffGene = arma::find(adjPTTest < ALPHA_THRESH * 2);
            corrP1 = corr(TF.cols(0, n1-1).t(), TG1.t());
            corrP2 = corr(TF.cols(n1, n1+n2-1).t(), TG2.t());
            pCombine = arma::min(corrP1, corrP2);
            tempIdx = arma::find(pCombine < ALPHA_THRESH);
            netIdx = arma::umat(tempIdx.n_elem, 2, arma::fill::zeros);
            convertIdxToRowCol(tempIdx, netIdx, pCombine.n_rows);

            for (unsigned int j = 0; j < diffGene.n_elem; ++j) {
                OTF1 = BO1.col(diffGene.at(j)) % TF.cols(0, n1 - 1);
                OTF2 = BO2.col(diffGene.at(j)) % TF.cols(n1, n1+n2-1);
                temp = netIdx.rows(arma::find(netIdx.col(1) == diffGene.at(j)));
                id1 = temp.col(0);
                id = arma::find(arma::sum(arma::abs(OTF1.t())) + arma::sum(arma::abs(OTF2.t())) > 0);
                id = arma::intersect(id1, id);

                for (unsigned int i = 0; i < id.n_elem; ++i) {
                    X1 = arma::join_horiz(OTF1.row(id.at(i)).t(), TG1.row(diffGene.at(j)).t());
                    X2 = arma::join_horiz(OTF2.row(id.at(i)).t(), TG2.row(diffGene.at(j)).t());
                    bivarRetTup = bivariateNormalConditionalLR(X1, X2);
                    LRi = {(double) id.at(i), (double) diffGene.at(j), std::get<3>(bivarRetTup), std::get<2>(bivarRetTup)};
                    LRSummary = arma::join_vert(LRSummary, LRi);
                }
            }
            // remove the first filler row added previously
            LRSummary.shed_row(0);
            // find NaN and +-Inf
            LRSummary.shed_rows(arma::find_nonfinite(LRSummary.col(2)));
            // m.elem((colInd - 1) * m.n_rows + (rowInd - 1));
            LRSummary.elem(LRSummary.n_rows + (arma::find(LRSummary.col(2) < tenSixteenthPow) - 1)).fill(tenSixteenthPow);

            LRSColTwo = LRSummary.col(2);
            phat = gammaQuantileMatch(LRSColTwo, arma::regspace(0.1, 0.1, 0.2));
            boost::math::gamma_distribution<> gammaDist(phat.at(0), phat.at(1));
            pGamma = arma::vec(LRSummary.n_rows, arma::fill::zeros);
            for (unsigned int k = 0; k < LRSummary.n_rows; ++k) {
                pGamma.at(k) = 1 - boost::math::cdf(gammaDist, LRSColTwo.at(k));
            }
            adjPGamma = fdrBH(pGamma);
            LRSummary = arma::join_horiz(LRSummary, pGamma, adjPGamma);
            id = arma::find(adjPGamma < 0.1);
            id1 = arma::uvec(id.n_elem).fill(2);

            sortedIdx = arma::sort_index(LRSummary.elem(arr2ind(id, id1, LRSummary.n_rows)), "descend");
            id1 = id.elem(sortedIdx);
            id = arma::uvec(id1.n_elem, arma::fill::zeros);

            extractElems(TFName, arma::conv_to<arma::uvec>::from(LRSummary.elem(arr2ind(id1, id, LRSummary.n_rows))));
            id.ones();
            extractElems(symbol, arma::conv_to<arma::uvec>::from(LRSummary.elem(arr2ind(id1, id, LRSummary.n_rows))));
            arma::mat LRTempMat = LRSummary.cols(2, 5);
            LRTempMat = LRTempMat.rows(id1);
        }
    } catch (...) {
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
        arma::umat idxMat = arma::umat(rrPos.n_elem, 2, arma::fill::zeros);
        convertIdxToRowCol(rrPos, idxMat, rrNRows);
        //row
        a.col(0) = arma::conv_to<arma::vec>::from(idxMat.col(0));
        // col
        a.col(1) = arma::conv_to<arma::vec>::from(idxMat.col(1));
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
        // TODO: sort will cause problem!!!
//        std::sort(symbol1.begin(), symbol1.end());
        std::vector<std::string>::iterator s1Start= symbol1.begin();
//        std::sort(symbol2.begin(), symbol2.end());
        std::vector<std::string>::iterator s2Start= symbol2.begin();

        arma::uvec f2 = arma::uvec(symbol.size(), arma::fill::zeros);
        arma::uvec f1 = arma::uvec(symbol.size(), arma::fill::zeros);
        unsigned int ind = 0;
        for (auto & t : symbol) {
            // using lower_bound because t is always found in symbol1 and symbol2
            f1(ind) = std::lower_bound(symbol1.begin(), symbol1.end(), t) - s1Start;
            f2(ind) = std::lower_bound(symbol2.begin(), symbol2.end(), t) - s2Start;
            ++ind;
        }

        arma::mat E1Mean = arma::mat(f1.n_rows, K1, arma::fill::zeros);
        arma::mat E2Mean = arma::mat(f2.n_rows, K2, arma::fill::zeros);

        double avg = 0.;
        unsigned int f1Idx;
        for (unsigned int i = 0; i < K1; ++i) {
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
        for (unsigned int i = 0; i < K2; ++i) {
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

//        std::sort(peakName1.begin(), peakName1.end());
        // std::vector<std::string>::iterator p1Start= peakName1.begin();
//        std::sort(peakName2.begin(), peakName2.end());
        // std::vector<std::string>::iterator p2Start= peakName2.begin();

        f2 = arma::uvec(peakNameIntersect2.size(), arma::fill::zeros);
        f1 = arma::uvec(peakNameIntersect1.size(), arma::fill::zeros);
        ind = 0;
        for (std::vector<std::string>::iterator t=peakNameIntersect1.begin(); t!=peakNameIntersect1.end(); ++t) {
            f1(ind) = std::lower_bound(peakName1.begin(), peakName1.end(), *t) - s1Start;
            ++ind;
        }
        ind = 0;
        for (std::vector<std::string>::iterator t=peakNameIntersect2.begin(); t!=peakNameIntersect2.end(); ++t) {
            f2(ind) = std::lower_bound(peakName2.begin(), peakName2.end(), *t) - s2Start;
            ++ind;
        }

        // need to sort here for binary search, but won't affect index since peakNameIntersect1 is not used later
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
        for (unsigned int i = 0; i < K1; ++i) {
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
        for (unsigned int i = 0; i < K2; ++i) {
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

