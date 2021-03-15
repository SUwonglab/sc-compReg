\name{cnmf_tsne}
\alias{cnmf_tsne}
\title{Joint clustering of scRNA-seq and scATAC-seq via Coupled Nonnegative Matrix Factorization}
\usage{
\method{cnmf_tsne}{default}(H1, H2, save.plot=T, path='./', ...)}
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
    \item{H1}{\code{matrix}. Output matrix from \code{cnmf}}
    \item{H2}{\code{matrix}. Output matrix from \code{cnmf}}
    \item{save.plot}{\code{logical} (default = T). Logical indicating whether or not the generated plots should be saved.}
    \item{path}{\code{character}. Path to the directory where the generated plots should be saved.}
}

\examples{
    cnmf_tsne(H1, H2, save.plot=T, path='./')
}
