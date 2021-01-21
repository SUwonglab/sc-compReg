library(scCompReg)
library(R.matlab)
library(Matrix)
library(tictoc)

pni = comp_reg_preprocess('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', token='\t')
s2 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample2.mat')
s1 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample1.mat')

tic()
output = cluster.profile(Matrix(s1$O1, sparse=T),
               Matrix(s1$E1, sparse=T),
                as.integer(s1$O1.idx) - 1,
               as.integer(s1$E1.idx) - 1,
               as.vector(unlist(s1$Symbol1, use.names = F)),
                as.vector(unlist(s1$PeakName1, use.names = F)),
               Matrix(s2$O2, sparse=T),
               Matrix(s2$E2, sparse=T),
               as.integer(s2$O2.idx) - 1,
               as.integer(s2$E2.idx) - 1,
               as.vector(unlist(s2$Symbol2, use.names = F)),
               as.vector(unlist(s2$PeakName2, use.names = F)),
                as.vector(pni$vo),
               as.vector(pni$vt)
               )
toc()

# K1 = max(max(s1$O1.idx), max(s1$E1.idx))
# K2 = max(max(s2$O2.idx), max(s2$E2.idx))
# symbol = intersect(s1$Symbol1, s2$Symbol2)
# f1 = match(symbol, s1$Symbol1)
# f2 = match(symbol, s2$Symbol2)
# E1.mean = matrix(NA, length(f1), K1)
# for (i in 1:K1) {
#     E1.mean[, i] = apply(s1$E1[f1, s1$E1.idx == i], 1, mean)
# }

write.csv(E1.mean, 'E1mean_r.csv')


new.output = subpopulation.link(output$E1.mean,
                                output$E2.mean,
                                output$O1.mean,
                                output$O2.mean)
write.csv(output$E1.mean, 'E1mean.csv')
write.csv(output$E2.mean, 'E2mean.csv')
write.csv(output$O1.mean, 'O1mean.csv')
write.csv(output$O2.mean, 'O2mean.csv')
