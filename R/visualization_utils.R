cnmf.tsne <- function(H1, H2,
                      save.plot=T,
                      path='./') {
    S10 <- apply(H1, 2, function(x) {return(which(x == max(x))[1])})
    S20 <- apply(H2, 2, function(x) {return(which(x == max(x))[1])})
    clust.assignments <- cbind(S10, S20)
    HH <- cbind(S10, S20)
    HH <- sweep(HH, 2, sqrt(rowSums(HH^2)), '/')
    HH.cos.sim <- HH / sqrt(rowSums(HH * HH))
    HH.cos.sim <- HH.cos.sim %*% t(HH.cos.sim)
    tsne.output <- Rtsne(HH.cos.sim, perplexity = 100,
                         check_duplicates = F,
                         is_distance = T,
                         pca_scale = T,
                         eta = 250)

    par(mfrow=c(2,1))
    num.unique.clust <- length(unique(clust.assignments))
    color.palette <- hue_pal()(num.unique.clust)
    color.assignment = clust.assignments
    for (j in 1:num.unique.clust) {
        color.assignment[num.unique.clust==j] = color.palette[j]
    }

    plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
         col=color.assignment, pch=19,
         cex=0.8, main = 't-SNE plot for RNA-seq and ATAC-seq')

    plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
         col=c(rep('red', length(S10)),
         rep('blue', length(S20))), pch=19,
         cex=0.8, main = 't-SNE plot for RNA-seq and ATAC-seq')
    legend("topleft", inset=.05,
           c("RNA-seq", "ATAC-seq"), fill=c('red', 'blue'), horiz=F, cex=0.7)

    if (save.plot) {
        par(mfrow=c(1,1))
        png(paste(path, 'joint_clusters_plot.png', sep='_'), units="in", width=8, height=8, res=300)
        plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
             col=color.assignment, pch=19,
             cex=0.8, main = 't-SNE plot for RNA-seq and ATAC-seq')
        dev.off()

        png(paste(path, 'joint_data_type_plot.png', sep='_'), units="in", width=8, height=8, res=300)
        plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
             col=color.assignment, pch=19,
             cex=0.8, main = 't-SNE plot for RNA-seq and ATAC-seq')
        legend("topleft", inset=.05,
               c("RNA-seq", "ATAC-seq"), fill=c('red', 'blue'), horiz=F, cex=0.7)
        dev.off()
    }
}
