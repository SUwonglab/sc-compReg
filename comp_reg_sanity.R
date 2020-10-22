pni = comp_reg_preprocess('/Users/Sophia/Desktop/BioStats/compreg/PeakName_intersect.txt', token='\t')
library(R.matlab)
library(Matrix)
s1 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample1.mat')
s2 = readMat('/Users/Sophia/Desktop/BioStats/compreg/sample2.mat')
clusterProfile(Matrix(s1$O1, sparse=T),
               Matrix(s1$E1, sparse=T),
                as.integer(s1$O1.idx) - 1,
               as.integer(s1$E1.idx) - 1,
               unlist(s1$Symbol1, use.names = F),
                unlist(s1$PeakName1, use.names = F),
               Matrix(s2$O2, sparse=T),
               Matrix(s2$E2, sparse=T),
               as.integer(s2$O2.idx) - 1,
               as.integer(s2$E2.idx) - 1,
               unlist(s2$Symbol2, use.names = F),
               unlist(s2$PeakName2, use.names = F),
                pni$vo,
               pni$vt
               )

library(scCompReg)

as.integer(s1$O1.idx)








