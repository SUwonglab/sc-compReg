library(scCompReg)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path = './example_data/'

cnmf.data <- readRDS(paste(path, 'cnmf_data.rds', sep=''))
cnmf.output <- cnmf(cnmf.data$PeakO,
                    cnmf.data$X,
                    cnmf.data$D,
                    k=3,
                    alpha=0.5,
                    beta_max_scale=5,
                    verbose=T)

cnmf_tsne(cnmf.output$H1, cnmf.output$H2, save.plot=F)

if (!require("Matrix")) install.packages("Matrix")
library(Matrix)
# log-transforming E matrix and O matrix
# sc-compReg requires log-transformed matrices
cnmf.data$PeakO@x = log2(cnmf.data$PeakO@x + 1)
writeMM(cnmf.data$PeakO,
        file=paste(path, 'O.mtx', sep=''))
cnmf.data$X@x = log2(cnmf.data$X@x + 1)
writeMM(cnmf.data$X,
        file=paste(path, 'E.mtx', sep=''))

write.table(cnmf.output$atac_cluster,
            paste(path, 'atac_cluster.txt', sep=''),
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t')

write.table(cnmf.output$rna_cluster,
            paste(path, 'rna_cluster.txt', sep=''),
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t')



