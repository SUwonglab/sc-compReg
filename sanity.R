library(scCompReg)
library(tictoc)
library(R.matlab)


file.mat <- readMat( '/Users/Sophia/Desktop/BioStats/RAd4_sc_data.mat')
X <- Matrix(file.mat[['X']], sparse=TRUE)
A <- Matrix(file.mat[['D']], sparse=TRUE)
peakO <- Matrix(log10(file.mat[['PeakO']] + 1), sparse=TRUE)

tic('cnmf')
output = cnmf(peakO,
         X,
         A,
         k=3,
         alpha=0.5,
         beta_max_scale=5)
toc()



e.symbol <- unlist(file.mat[[mat_obj$e.sym]], use.names=FALSE)
p.symbol <- unlist(file.mat[[mat_obj$p.sym]], use.names=FALSE)


output = cnmf(mat_path = '/Users/Sophia/Desktop/BioStats/RAd4_sc_data.mat',
              k_components=3,
              mat_obj=list('A'='D',
                           'po'='PeakO',
                           'X'='X',
                           'e.sym'='Symbol',
                           'p.sym'='PeakName'),
              lambda_1 = 0.1688,
              lambda_2= 0.0018,
              mu=500)
toc()

length(which(s2w[1:464] == s2[1:464]))
length(which(s1w[1:415] == s1))

length(which(s1w[1:415] == j))

length(which(j == s1))
hist(s1w[1:415])

j = output$S1
j[j==1]=4
j[j==3]=1
j[j==4]=3
hist(j)

hist(s1)
hist(s2)
hist(j2)
j2 = output$S2
j2[j2==2]=4
j2[j2==1]=2
j2[j2==4]=1

length(which(j2==s2))
hist(j2)
hist(s2)


order(c(length(s2[s2==1]), length(s2[s2==2]), length(s2[s2==3])))
hist(j2)
j2[j2==2]=4
j2[j2]



length(which(j2==s2))
length(which(j2==s2w[1:464]))
length(which(s2w[1:464] == s2))


length(s1)





system.time({output <- cnmf(peak_O_path = "/Users/Sophia/Desktop/BioStats/PeakO.txt",
              p_symbol_path = '/Users/Sophia/Desktop/BioStats/PeakName.txt',
              e_symbol_path = "/Users/Sophia/Desktop/BioStats/symbol.txt",
              e_matrix_path = "/Users/Sophia/Desktop/BioStats/E.txt",
              pe_path = "/Users/Sophia/Desktop/BioStats/peak_gene_100k_corr.bed",
              k_components = 2,
              lambda_2=10, o.loss='mse')})


install.packages('M3C')
library(M3C)
umap(output$)

























library(uwot)
m = umap(t(output$H1), ret_model = TRUE)
library(ggplot2)
df = data.frame(m$embedding)
df$cluster = as.factor(s1w[1:200])
ggplot(df, aes(x=X1, y=X2)) + geom_point(aes(colour=cluster), size=3) + labs(y='UMAP_2', x='UMAP_1', fill='Cell Cluster') + theme(legend.title = element_text(color = "black", size = 10),
            legend.text = element_text(color = "darkgrey"), plot.title=element_text(face="bold",hjust=0.5)) + theme_minimal() + ggtitle('UMAP for H_1 Matrix')
ggsave('./POH.png', plot = last_plot())

m = umap(t(output$H2), ret_model = TRUE)
library(ggplot2)
df = data.frame(m$embedding)
df$cluster = as.factor(s2w[1:200])
ggplot(df, aes(x=X1, y=X2)) + geom_point(aes(colour=cluster), size=3) + labs(y='UMAP_2', x='UMAP_1', fill='Cell Cluster') + theme(legend.title = element_text(color = "black", size = 10),
                                                                                                                                  legend.text = element_text(color = "darkgrey"), plot.title=element_text(face="bold",hjust=0.5)) + theme_minimal() + ggtitle('UMAP for H_2 Matrix')
ggsave('./EH.png', plot = last_plot())

install.packages('Seurat')

output = cnmf(peak_O_path = "/Users/Sophia/Desktop/BioStats/PeakO.txt",
              p_symbol_path = '/Users/Sophia/Desktop/BioStats/PeakName.txt',
              e_symbol_path = "/Users/Sophia/Desktop/BioStats/symbol.txt",
              e_matrix_path = "/Users/Sophia/Desktop/BioStats/E.txt",
              pe_path = "/Users/Sophia/Desktop/BioStats/peak_gene_100k_corr.bed",
              k_components = 2,
              lambda_2=30)

output$S1

hist(s1)
hist(output$S1)
hist(s2)
hist(output$S2)
s1 = as.matrix(read.table('/Users/Sophia/Desktop/BioStats/Label_ATACseq', sep='\t', header=FALSE))
s2 = as.matrix(read.table("/Users/Sophia/Desktop/BioStats/Label_RNAseq", sep='\t', header=FALSE))

library(R.matlab)
test=readMat('/Users/Sophia/Desktop/BioStats/RAd4_sc_data.mat')

PO = as.matrix(read.table("/Users/Sophia/Desktop/BioStats/PeakO.txt", header=FALSE))
E = as.matrix(read.table("/Users/Sophia/Desktop/BioStats/E.txt", header=FALSE))

## MINE
os1 = t(output$S1)
os2 = t(output$S2)
var(as.vector(PO[, which(os1==1)])) + var(as.vector(PO[, which(os1==2)]))
var(as.vector(E[, which(os2==1)])) + var(as.vector(E[, which(os2==2)]))

## Wanwen's
s1w = as.matrix(read.table('/Users/Sophia/Desktop/BioStats/CoupledNMF/scATAC-result.txt', sep='\t', header=FALSE))
s2w = as.matrix(read.table("/Users/Sophia/Desktop/BioStats/CoupledNMF/scRNA-result.txt", sep='\t', header=FALSE))
hist(s1)
length(s1)
var(as.vector(PO[, which(s1==1)])) + var(as.vector(PO[, which(s1==2)]))
var(as.vector(E[, which(s2==1)])) + var(as.vector(E[, which(s2==2)]))

# ## Random baseline
# if (!require('extraDistr')) {install.packages('extraDistr')}
# library(extraDistr)
# r1 =rdunif(200, 1, 2)
# r2 = rdunif(200, 1, 2)
# var(as.vector(PO[, which(r1==1)])) + var(as.vector(PO[, which(r1==2)]))
# var(as.vector(E[, which(r2==1)])) + var(as.vector(E[, which(r2==2)]))
output$S2

sum(which(os1==1)>100) + sum(which(os1==2) <100)
sum(which(os2==1)<100) + sum(which(os2==2) >100)

sum(which(s1w==1)<100) + sum(which(s1w==2) >100)
sum(which(s2w==1)<100) + sum(which(s2w==2) >100)

write.table(output$S1,file="S1.csv", row.names=FALSE,sep=',', col.names=FALSE)
write.table(output$S2,file="S2.csv", row.names=FALSE,sep=',', col.names=FALSE)


write.table(as.matrix(test$D),file='D.csv',row.names = FALSE,sep=',',col.names = FALSE)
write.table(test$PeakO,file='PO.csv',row.names = FALSE,sep=',',col.names = FALSE)
write.table(test$X,file='X.csv',row.names = FALSE,sep=',',col.names = FALSE)
write.table(noquote(unlist(test$PeakName, use.names=FALSE)),file='PN.txt',row.names = FALSE,sep=' ',col.names = FALSE, quote = FALSE)
write.table(noquote(unlist(test$Symbol, use.names=FALSE)),file='symb.txt',row.names = FALSE,sep=' ', col.names = FALSE, quote=FALSE)




length(which(s2 != os2))
length(which(s1 != os1))



