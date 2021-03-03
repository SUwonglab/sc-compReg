library(scCompReg)
library(Matrix)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path = './data/'

cnmf.data <- readRDS(paste(path, 'cnmf_data.rds', sep=''))
cnmf.output <- cnmf(cnmf.data$PeakO,
                    cnmf.data$X,
                    cnmf.data$D,
                    k=3,
                    alpha=0.5,
                    beta_max_scale=5,
                    verbose=T)

cnmf_tsne(cnmf.output$H1, cnmf.output$H2, save.plot=F)

write.table(cnmf.output$atac_cluster,
            paste(path, 'atac_cluster.txt', sep=''),
            row.names = F,
            quote = F,
            sep = '\t')

write.table(cnmf.output$rna_cluster,
            paste(path, 'atac_cluster.txt', sep=''),
            row.names = F,
            quote = F,
            sep = '\t')

