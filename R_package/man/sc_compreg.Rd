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
    \item{O1}{\code{matrix} or \code{sparseMatrix}. Peak by cell matrix of scATAC-seq from sample 1 (healthy).}
    \item{E1}{\code{matrix} or \code{sparseMatrix}. Gene by cell matrix of scRNA-seq from sample 1 (healthy).}
    \item{O1.idx}{\code{numeric} or \code{matrix} or \code{integer}. Clustering labels of the \code{O1} matrix columns from sample 1.}
    \item{E1.idx}{\code{numeric} or \code{matrix} or \code{integer}. Clustering labels of the \code{E1} matrix columns from sample 1.}
    \item{symbol1}{\code{character}. symbol1 is a vector of characters, indicating the raw names of E1 from sample 1.}
    \item{peak.name1}{\code{character} peak.name1 is a vector of characters, indicating the raw names of O1 from sample 1.}
    \item{O2}{\code{matrix} or \code{sparseMatrix}. Peak by cell matrix of scATAC-seq from sample 2 (diseased).}
    \item{E2}{\code{matrix} or \code{sparseMatrix}. Gene by cell matrix of scRNA-seq from sample 1 (diseased).}
    \item{O2.idx}{\code{numeric} or \code{matrix} or \code{integer}. Clustering labels of the \code{O2} matrix columns from sample 2.}
    \item{E2.idx}{\code{numeric} or \code{matrix} or \code{integer}. Clustering labels of the \code{E2} matrix columns from sample 2.}
    \item{symbol2}{\code{character}. symbol2 is a vector of characters, indicating the raw names of E2 from sample 2.}
    \item{peak.name2}{\code{character} peak.name2 is a vector of characters, indicating the raw names of O2 from sample 2.}
    \item{motif.name}{\code{list} or \code{character}. A vector of characters indicating the name of motifs.}
    \item{motif.weight}{\code{numeric} or \code{matrix}. motif.weight is a vector of numeric indicating the weight of motif based on the number of their matching sites on the genome.}
    \item{motif.file}{\code{list} motif.file is the loaded motif_target file (\code{list}) returned from \code{mfbs_load}.}
    \item{peak.name.intersect.dir}{\code{character}. Path to the directory of the peak_name_intersect file. peak_name_intersect file includes overlapped peaks from two samples.}
    \item{peak.gene.prior.dir}{\code{character}. Path to the directory of the peak_gene_prior file. peak_gene_prior file includes prior of peak-gene association from PECA model, which is used to constrcut the coupling matrix D.}
    \item{sep.char}{\code{character}. The delimiter or separating character in the file, peak_name_intersect. Default: \code{'\\t'}}
    \item{thresh}{\code{numeric}.}
    \item{sig.level}{\code{numeric}.}
    \item{num.top.tf}{\code{numeric}.}
    \item{d0.default}{\code{numeric}.}
    \item{verbose}{\code{logical}.}
}

\value{
    \item{n.pops}{\code{numeric}. Number of linked subpopulations (cluster matching) between sample 1 and sample 2.}
    \item{match}{\code{matrix} Subpopulation matching table between sample 1 and sample 2.}
    \item{diff.net}{\code{list}. Differential gene regulatory network between two linked subpopulation.}
    \item{hub.tf}{\code{list}. TF names and their in-degree in each of the differentional regulatory network.}
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
