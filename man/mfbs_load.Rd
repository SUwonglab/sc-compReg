\name{mfbs_load}
\alias{mfbs_load}
\title{Single cell comparative regulatory analysis based on scRNA-seq and scATAC-seq from two different conditions}
\usage{
\method{mfbs_load}{default}(motif.target.dir, ...)
}
\description{
Efficiently loads the \code{motif_target} file necessary for the \code{sc_compreg} methods via \code{c++}.}
\author{Sophia Lu \email{sophialu@stanford.edu}, Zhana Duren}

\details{
For details of the implementation and method, please see Duren et al. (2020)
}

\references{
\emph{Comparative regulatory analysis of single cell data reveals a novel B cell subpopulation in chronic lymphocytic leukemia}. (2020)
Zhana Duren, Wenhui Sophia Lu, Joseph G. Arthur, Preyas Shah, Jingxue Xin,  Francesca Meschi, Miranda Lin Li, Corey M. Nemec, Yifeng Yin, and Wing Hung Wong
}

\arguments{
    \item{motif.target.dir}{\code{character}. Path to the directory where \code{motif_target} file is stored.}
}

\value{
    \item{motif.file}{\code{list}. A loaded list to be passed to \code{sc_compreg}.}
}


\examples{
    motif.file = mfbs_load(paste(path, motif.target.dir, sep=''))
}
