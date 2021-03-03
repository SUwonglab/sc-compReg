\name{cnmf}
\alias{cnmf}
\title{Joint clustering of scRNA-seq and scATAC-seq via Coupled Nonnegative Matrix Factorization}
\usage{
\method{cnmf}{default}(peak.o, X, D, k, beta.min, alpha=0.5, beta.max.scale=5, verbose=T, ...)
}
\description{
Performs joint clustering of scRNA-seq and scATAC-seq via coupled nonnegative matrix factorization.}
\author{Sophia Lu \email{sophialu@stanford.edu}, Zhana Duren}

\details{
For details of the implementation and method, please see Duren et al. (2020)
}

\references{
\emph{Comparative regulatory analysis of single cell data reveals a novel B cell subpopulation in chronic lymphocytic leukemia}. (2020)
Zhana Duren, Wenhui Sophia Lu, Joseph G. Arthur, Preyas Shah, Jingxue Xin,  Francesca Meschi, Miranda Lin Li, Corey M. Nemec, Yifeng Yin, and Wing Hung Wong
}

\arguments{
    \item{peak.o}{\code{matrix} or \code{sparseMatrix}. Matrix of the openness counts (ATAC-seq data).}
    \item{X}{\code{matrix} or \code{sparseMatrix}. Matrix of gene expression (RNA-seq data).}
    \item{D}{\code{matrix} or \code{sparseMatrix}.Coupling matrix learned from paired expression and chromatin accessibility data across diverse panel of cellular contexts.}
    \item{k}{\code{numeric} or \code{integer}. Number of clusters to identify.}
    \item{beta.min}{\code{numeric}. If \code{beta.min} is missing, then it is set to
    \code{10.0^seq(beta.max.scale, -3, -1)}.}
    \item{alpha}{\code{numeric} Tuning parameter balancing the weight of scRNA-seq and scATAC-seq. Default: 0.5.}
    \item{beta.max.scale}{\code{numeric}.}
}

\value{
    \item{W1}{\code{matrix}. Final W1 matrix.}
    \item{W2}{\code{matrix}. Final W2 matrix.}
    \item{H1}{\code{matrix}. Final H1 matrix.}
    \item{H2}{\code{matrix}. Final H2 matrix.}
    \item{lambda1}{\code{numeric}.}
    \item{lambda2}{\code{numeric}.}
    \item{score}{\code{numeric}.}
    \item{cluster1}{\code{vector}. A vector indicating cluster assignments for ATAC-seq data.}
    \item{cluster2}{code{vector}. A vector indicating cluster assignments for RNA-seq data.}
}


\examples{
    library(scCompReg)
    library(R.matlab)

    file.mat <- readMat( '~/data/cnmf_data.mat')
    X <- Matrix(file.mat[['X']], sparse=TRUE)
    D <- Matrix(file.mat[['D']], sparse=TRUE)
    peak.o <- Matrix(log10(file.mat[['PeakO']] + 1), sparse=TRUE)

    output = cnmf(peak.o,
                  X,
                  D,
                  k=3,
                  alpha=0.5,
                  beta.max.scale=5,
                  verbose=T)

}
