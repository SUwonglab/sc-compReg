library(scCompReg)

rm(list = ls())

# set the path to the document with all of the input files
# path should end in '/'
# for example, something like path = 'Desktop/compreg_input/'
path <- './example_data/'

# # to convert the previously generated ```peak_gene_coupling_matrix.txt``` by
# # ```cnmf_process_data.sh``` to D matrix
# # change the directory to where ```peak_gene_coupling_matrix.txt```
# # is stored and run the lines below
# preprocess.path <- './preprocess_data/'
# # load peak.name
# # load symbol
# symbol <- as.vector(unlist(read.table(paste(path, 'cnmf_symbol.txt', sep=''),
#                                          header=F), use.names=F))
# peak.name <- as.vector(unlist(read.table(paste(path, 'cnmf_peak_name.txt', sep=''),
#                                       header=F), use.names=F))
# print('Printing first 5 rows of sample 1 peak name...')
# print(symbol[1:5])
# print('Printing first 5 rows of sample 1 symbol...')
# print(peak.name[1:5])
# D <- cnmf_load_coupling_matrix(paste(preprocess.path, 'peak_gene_coupling_matrix.txt', sep=''),
#                                peak.name,
#                                symbol)


cnmf.data <- readRDS(paste(path, 'cnmf_data.rds', sep=''))
cnmf.output <- cnmf(cnmf.data$PeakO,
                    cnmf.data$X,
                    cnmf.data$D,
                    k=3,
                    alpha=0.5,
                    beta_max_scale=5,
                    verbose=F)

cnmf_tsne(cnmf.output$H1, cnmf.output$H2, path=path, save.plot=F, perplexity=100)

# if data is not log-transformed, perform log-transformation
# if (!require("Matrix")) install.packages("Matrix")
# library(Matrix)
# # log-transforming E matrix and O matrix
# # ```cnmf```` requires log-transformed matrices
# cnmf.data$PeakO@x <- log2(cnmf.data$PeakO@x + 1)
# writeMM(cnmf.data$PeakO,
#         file=paste(path, 'O.mtx', sep=''))
# cnmf.data$X@x <- log2(cnmf.data$X@x + 1)
# writeMM(cnmf.data$X,
#         file=paste(path, 'E.mtx', sep=''))

write.table(cnmf.output$atac_cluster,
            paste(path, 'atac_cluster.txt', sep=''),
            row.names = T,
            col.names = F,
            quote = F,
            sep = '\t')

write.table(cnmf.output$rna_cluster,
            paste(path, 'rna_cluster.txt', sep=''),
            row.names = T,
            col.names = F,
            quote = F,
            sep = '\t')



