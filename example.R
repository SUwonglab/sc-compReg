library(scCompReg)
library(R.matlab)
library(Matrix)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path = ''
s1 = readMat(paste(path, 'sample1.mat', sep=''))
s2 = readMat(paste(path, 'sample2.mat', sep=''))
O1 = s1$O1
O1.idx = s1$O1.idx
symbol1 = s1$Symbol1
O2 = s2$O2
O2.idx = s2$O2.idx
symbol2 = s2$Symbol2
peak.name1 = s1$PeakName1
E1.idx = s1$E1.idx
E1 = s1$E1
E2 = s2$E2
E2.idx = s2$E2.idx
peak.name2 = s2$PeakName2
peak.name.intersect.dir = paste(path, 'PeakName_intersect.txt', sep='')
motif.target.dir = paste(path, 'MotifMatch_human_rmdup.mat', sep='')
motif.target.dir = paste(path, 'MotifTarget.txt', sep='')
peak.gene.prior.dir = paste(path, 'peak_gene_prior_intersect.bed', sep='')
sep.char=' '

motif.mat = readMat(motif.mat.dir)
motif.name = unlist(motif.mat$motifName, use.names=F)
motif.weight = as.numeric(unlist(motif.mat$motifWeight, use.names=F))
match2 = unlist(motif.mat$Match2, use.names = F)
m2.half.idx = length(match2) / 2
match2 = list('a' = match2[1: m2.half.idx],
              'b' = match2[(m2.half.idx + 1):length(match2)])

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
    write.table(output$hub.tf[[i]], paste(path, 'tf_', i, '.txt', sep=''),
                col.names = F)
    write.table(output$diff.net[[i]], paste(path, 'diff_net_', i, '.txt', sep=''),
                col.names = F)
}

