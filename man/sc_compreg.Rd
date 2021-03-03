\name{sc_compreg}
\alias{sc_compreg}
\title{Single cell comparative regulatory analysis based on scRNA-seq and scATAC-seq from two different conditions}
\usage{
\method{sc_compreg}{default}(O1, E1, O1.idx, E1.idx, symbol1, peak.name1, O2, E2, O2.idx, E2.idx, symbol2, peak.name2, motif.name, motif.weight, match2, motif.file, peak.name.intersect.dir, peak.gene.prior.dir, sep.char = '\t', thresh = 0.2, sig.level = 0.05, num.top.tf = 5000L, d0.default = 500000L, verbose = TRUE, ...)
}
\description{
\code{sc_compreg} takes scRNA-seq and scATAC-seq data matrix from two samples (sample 1 from healthy individual, sample 2 from diagnosed individual) as well as their clusteirng labels as inputs, performs subpopulation link across two sample, and detects differential regulatory analysis between each pair of linked subpopulations.}

\author{Sophia Lu \email{sophialu@stanford.edu}, Zhana Duren}

\note{\code{sample 1} is the sample from a healthy individual, and \code{sample 2} is the sample from a diagnosed individual.}

\details{
For details of the implementation and method, please see Duren et al. (2020)
}

\references{
\emph{Comparative regulatory analysis of single cell data reveals a novel B cell subpopulation in chronic lymphocytic leukemia}. (2020)
Zhana Duren, Wenhui Sophia Lu, Joseph G. Arthur, Preyas Shah, Jingxue Xin,  Francesca Meschi, Miranda Lin Li, Corey M. Nemec, Yifeng Yin, and Wing Hung Wong
}

\arguments{
    \item{O1}{\code{matrix} or \code{sparseMatrix}. Matrix of the openness counts (ATAC-seq data).}
    \item{E1}{\code{matrix} or \code{sparseMatrix}. Matrix of the gene expression (RNA-seq data).}
    \item{O1.idx}{\code{matrix} or \code{sparseMatrix}. }
    \item{E1.idx}{\code{matrix} or \code{sparseMatrix}.}
    \item{symbol1}{\code{matrix} or \code{sparseMatrix}.}
    \item{peak.name1}{\code{matrix} or \code{sparseMatrix}.}
    \item{O2}{\code{matrix} or \code{sparseMatrix}.}
    \item{E2}{\code{matrix} or \code{sparseMatrix}.}
    \item{O2.idx}{\code{matrix} or \code{sparseMatrix}.}
    \item{E2.idx}{\code{matrix} or \code{sparseMatrix}.}
    \item{symbol2}{\code{matrix} or \code{sparseMatrix}.}
    \item{peak.name2}{\code{matrix} or \code{sparseMatrix}.}
    \item{motif.name}{\code{matrix} or \code{sparseMatrix}.}
    \item{motif.weight}{\code{matrix} or \code{sparseMatrix}.}
    \item{motif.file}{\code{matrix} or \code{sparseMatrix}.}
    \item{peak.name.intersect.dir}{\code{matrix} or \code{sparseMatrix}.}
    \item{peak.gene.prior.dir}{\code{matrix} or \code{sparseMatrix}.}
    \item{sep.char}{\code{matrix} or \code{sparseMatrix}.}
    \item{thresh}{\code{matrix} or \code{sparseMatrix}.}
    \item{sig.level}{\code{matrix} or \code{sparseMatrix}.}
    \item{num.top.tf}{\code{matrix} or \code{sparseMatrix}.}
    \item{d0.default}{\code{matrix} or \code{sparseMatrix}.}
    \item{verbose}{\code{matrix} or \code{sparseMatrix}.}
}

\value{
    \item{W1}
}


\examples{
    library(scCompReg)

    path = './data/'
    sample1 = readRDS(paste(path, 'sample1.rds', sep = ''))
    sample2 = readRDS(paste(path, 'sample2.rds', sep = ''))

    peak.name.intersect.dir = paste(path, 'PeakName_intersect.txt', sep='')
    motif.target.dir = paste(path, 'MotifTarget.txt', sep='')
    peak.gene.prior.dir = paste(path, 'peak_gene_prior_intersect.bed', sep='')
    sep.char=' '

    motif = readRDS(paste(path, 'motif.rds', sep=''))
    motif.file = readRDS(paste(path, 'motif_file.rds', sep=''))

    compreg.output = sc_compreg(sample1$O1,
                                sample1$E1,
                                sample1$O1.idx,
                                sample1$E1.idx,
                                sample1$symbol1,
                                sample1$peak.name1,
                                sample2$O2,
                                sample2$E2,
                                sample2$O2.idx,
                                sample2$E2.idx,
                                sample2$symbol2,
                                sample2$peak.name2,
                                motif$motif.name,
                                motif$motif.weight,
                                motif$match2,
                                motif.file,
                                peak.name.intersect.dir,
                                peak.gene.prior.dir)
}
