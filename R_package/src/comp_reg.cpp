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
            Rcpp::checkUserInterrupt();
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


template <class T, class U>
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
std::tuple<arma::uvec, arma::vec> isMember(const T& vecA,
                                          const U& vecB) {
    auto BBegin = vecB.begin();
    auto BEnd = vecB.end();
    arma::uvec lia = arma::uvec(vecA.size(), arma::fill::zeros);
    arma::vec locb = arma::vec(vecA.size(), arma::fill::zeros);

    unsigned int insertIdx = 0;
    for (auto t : vecA) {
        auto p = std::find(BBegin, BEnd, t);
        if (p != BEnd) {
            lia.at(insertIdx) = 1;
            locb.at(insertIdx) = std::distance(BBegin, p);
        } else {
            lia.at(insertIdx) = 0;
            locb.at(insertIdx) = arma::datum::nan;
        }
        ++insertIdx;
    }
    return std::tuple<arma::uvec, arma::vec>{lia, locb};
}

template <class T, class U>
T intersection(T v1,
               U v2){
    T v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          std::back_inserter(v3));
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

// [[Rcpp::export]]
arma::vec threshK(const arma::mat& A,
                     unsigned int thresh) {
    try {
        arma::vec retVec = arma::vec(A.n_rows, arma::fill::zeros);
        arma::rowvec sortedVec;
        for (unsigned int i = 0; i < A.n_rows; ++i) {
            sortedVec = arma::sort(A.row(i), "descend");
            retVec.at(i) = sortedVec.at(thresh);
        }
        return retVec;
    } catch (...) {
        ::Rf_error("c++ exception");
    }
}

/** For fast rcpp sparse matrix multiplication.
 *
 * @param A - sp_mat
 * @param B - sp_mat
 * @return - A * B, sp_mat
 */
// [[Rcpp::export]]
arma::sp_mat mult(const arma::sp_mat& A,
                  const arma::sp_mat& B) {
    return A * B;
}


// [[Rcpp::export]]
arma::rowvec colMax(const arma::sp_mat& X) {
    arma::mat Y(X);
    return arma::max(Y, 0);
}


void parseMotifTarget(std::string filePath,
                      std::vector<std::string>& stringVec1,
                      std::vector<std::string>& stringVec2,
                      std::vector<float>& doubleVec,
                      char sep) {
    char buffer[BUFFER_SIZE];
    char *a = (char *)malloc(BUFFER_SIZE);
    char *b = (char *)malloc(BUFFER_SIZE);
    float c;
    int scanRet;
    FILE* f = fopen(filePath.c_str(), "r");
    std::string temp;
    std::string inputFormat = std::string("%s") +
                              sep +
                              std::string("%s") +
                              sep +
                              std::string("%f");
    const char *inputCStr = inputFormat.c_str();
    while (true) {
        if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
        scanRet = sscanf(buffer, inputCStr, a, b, &c);
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
Rcpp::List mfbsLoad(const std::string& motifTargetPath,
                    char sep) {
    try {
        std::vector<std::string> strVec1, strVec2;
        std::vector<float> floatVec3;
        unsigned int index;
        parseMotifTarget(motifTargetPath, strVec1, strVec2, floatVec3, sep);
        return List::create(Named("C1") = strVec1,
                            Named("C2") = strVec2,
                            Named("C3") = floatVec3);
    } catch (...) {
        ::Rf_error("c++ exception");
    }
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

// [[Rcpp::export]]
arma::mat corrTest(const arma::mat& X,
                   const arma::mat& Y) {
    try {
        arma::mat testStat = arma::cor(X, Y);
        double num = sqrt(X.n_rows - 2);
        testStat = testStat * num / arma::sqrt(1 - arma::pow(testStat, 2));
        return testStat;
    } catch(...) {
        ::Rf_error("c++ exception");
    }
}



// [[Rcpp::export]]
arma::vec accumArrayMin(const arma::vec& subs,
                        const arma::vec& val) {
    try{
        unsigned int maxSubs = (unsigned int) subs.max() + 1;
        unsigned int subsNRows = subs.n_rows;
        arma::vec ret = arma::vec(maxSubs, arma::fill::zeros);
        arma::uvec idx;
        for (unsigned int i = 0; i < maxSubs; ++i) {
            idx = arma::find(subs == i);
            if (idx.n_elem > 0) {
                ret.at(i) = arma::min(val.elem(idx));
            }
        }
        return ret;
    } catch(...) {
        ::Rf_error("c++ exception");
    }

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


arma::mat colwiseElemMult(arma::mat A, const arma::vec& B) {
    for (unsigned int i = 0; i < A.n_cols; ++i) {
        A.col(i) *= B.at(i);
    }
    return A;
}


// [[Rcpp::export]]
arma::sp_mat computeBZ(const arma::mat& tfBinding, const arma::vec& OM,
                    const arma::sp_mat& beta) {
    try {
        arma::mat rowMultTF = colwiseElemMult(tfBinding, OM);
        return arma::sp_mat(rowMultTF) * beta;
    } catch (...) {
        ::Rf_error("c++ exception");
    }
}


// [[Rcpp::export]]
Rcpp::List compRegLoad(const std::string& peakGenePriorPath) {
    try {
        std::vector<std::string> stringVec1, stringVec2;
        std::vector<float> floatVec3, floatVec4;
        char buffer[BUFFER_SIZE];
        char *a = (char *)malloc(BUFFER_SIZE);
        char *b = (char *)malloc(BUFFER_SIZE);
        float c, d;
        int scanRet;
        FILE* f = fopen(peakGenePriorPath.c_str(), "r");
        std::string temp;
        while (true) {
            if (fgets(buffer, BUFFER_SIZE, f) == NULL) break;
            scanRet = sscanf(buffer, "%s %s %f %f", a, b, &c, &d);
            temp = a;
            stringVec1.push_back(temp);
            temp = b;
            stringVec2.push_back(temp);
            floatVec3.push_back(c);
            floatVec4.push_back(d);
        }
        fclose(f);
        free(a);
        free(b);
        return Rcpp::List::create(Named("C1") = stringVec1,
                                  Named("C2") = stringVec2,
                                  Named("C3") = floatVec3,
                                  Named("C4") = floatVec4);
    } catch (...) {
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


