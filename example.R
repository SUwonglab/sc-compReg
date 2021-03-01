library(scCompReg)
library(R.matlab)
library(Matrix)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path = './data/'
sample1 = readRDS(paste(path, 'sample1.rds', ssep = ''))
sample2 = readRDS(paste(path, 'sample2.rds', ssep = ''))

peak.name.intersect.dir = paste(path, 'PeakName_intersect.txt', sep='')
motif.target.dir = paste(path, 'MotifTarget.txt', sep='')
peak.gene.prior.dir = paste(path, 'peak_gene_prior_intersect.bed', sep='')
sep.char=' '


motif = readRDS(paste(path, 'motif.rds', sep=''))

compreg.output = sc_compreg(O1,
                            E1,
                            O1.idx,
                            E1.idx,
                            symbol1,
                            peak.name1,
                            O2,
                            E2,
                            O2.idx,
                            E2.idx,
                            symbol2,
                            peak.name2,
                            motif.name,
                            motif.weight,
                            match2,
                            peak.name.intersect.dir,
                            motif.target.dir,
                            peak.gene.prior.dir)

for (i in 1:compreg.output$n.pops) {
    write.table(compreg.output$hub.tf[[i]], paste(path, 'tf_', i, '.txt', sep=''),
                row.names = F)
    write.table(compreg.output$diff.net[[i]], paste(path, 'diff_net_', i, '.txt', sep=''),
                row.names = F)
}

