#include "../inst/include/quantile_norm.h"


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
