library(scCompReg)
library(R.matlab)
library(Matrix)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path = './data/'
sample1 = readRDS(paste(path, 'sample1.rds', sep = ''))
sample2 = readRDS(paste(path, 'sample2.rds', sep = ''))

peak.name.intersect.dir = paste(path, 'PeakName_intersect.txt', sep='')
motif.target.dir = paste(path, 'MotifTarget.txt', sep='')
peak.gene.prior.dir = paste(path, 'peak_gene_prior_intersect.bed', sep='')
sep.char=' '

motif = readRDS(paste(path, 'motif.rds', sep=''))
# To load the MotifTarget file, use the following line
# motif.file = mfbs_load(paste(path, motif.target.dir, sep=''))

# To run through the example with example data
# simply run the line below
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

for (i in 1:compreg.output$n.pops) {
    write.table(compreg.output$hub.tf[[i]],
                paste(path, 'tf_', i, '.txt', sep=''),
                row.names = F,
                quote = F,
                sep = '\t')
    write.table(compreg.output$diff.net[[i]],
                paste(path, 'diff_net_', i, '.txt', sep=''),
                row.names = F,
                quote = F,
                sep = '\t')
}

