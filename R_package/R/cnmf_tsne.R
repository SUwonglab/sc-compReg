cnmf_tsne <- function(H1, H2,
                      save.plot=T,
                      path='./',
                      perplexity=100,
                      ...) {
    if (! is (H1, 'matrix')) {
        stop("H1 must be a matrix (output from cnmf).")
    }

    if (! is (H1, 'matrix')) {
        stop("H2 must be a matrix (output from cnmf).")
    }

    if (! is (save.plot, 'logical')) {
        stop("save.plot must be a logical indicating whether or not plots should be saved.")
    }

    if (! is (path, 'character')) {
        stop("path must be of type format indicating the directory where the plots should be saved.")
    }


    S10 <- apply(H1, 2, function(x) {return(which(x == max(x))[1])})
    S20 <- apply(H2, 2, function(x) {return(which(x == max(x))[1])})
    clust.assignments <- c(S10, S20)
    HH <- cbind(H1, H2)
    HH <- t(sweep(HH, 2, sqrt(colSums(HH^2)), '/'))
    HH.cos.sim <- HH / sqrt(rowSums(HH * HH))
    HH.cos.sim <- HH.cos.sim %*% t(HH.cos.sim)
    tsne.output <- Rtsne(HH.cos.sim, perplexity = perplexity,
                         check_duplicates = F,
                         is_distance = F,
                         pca_scale = T,
                         eta = 250)

    par(mfrow=c(1,2))
    num.unique.clust <- length(unique(clust.assignments))
    color.palette <- hue_pal()(num.unique.clust)
    color.assignment <- clust.assignments
    for (j in 1:num.unique.clust) {
        color.assignment[color.assignment==j] <- color.palette[j]
    }

    plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
         col=color.assignment, pch=19,
         cex=0.8, main = 't-SNE plot colored by clustering output',
         xlab='x', ylab='y')
    legend("topleft", inset=.05,
           title="Clusters",
           c(as.character(sort(unique(clust.assignments)))),
           fill=c(color.palette[sort(unique(clust.assignments))]),
           horiz=F, cex=0.7,
           bty = "n")

    plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
         col=c(rep('red', length(S10)),
         rep('blue', length(S20))), pch=19,
         cex=0.8, main = 't-SNE plot for RNA-seq and ATAC-seq',
         xlab='x', ylab='y')
    legend("topright", inset=.05,
           c("RNA-seq", "ATAC-seq"),
           fill=c('red', 'blue'), horiz=F, cex=0.7,
           bty = "n")

    if (save.plot) {
        par(mfrow=c(1,1))
        png(paste(path, 'joint_clusters_plot.png', sep=''), units="in", width=8, height=8, res=300)
        plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
             col=color.assignment, pch=19,
             cex=1, cex.main=2, cex.axis=2,
             main = 't-SNE plot colored by clustering output',
             xlab='x', ylab='y')
        legend("topright", inset=.05,
               title="Clusters",
               c(as.character(sort(unique(clust.assignments)))),
               fill=c(color.palette[sort(unique(clust.assignments))]),
               horiz=F, cex=2,
               bty = "n")
        dev.off()

        png(paste(path, 'joint_data_type_plot.png', sep=''), units="in", width=8, height=8, res=300)
        plot(tsne.output$Y[, 1], tsne.output$Y[, 2],
             col=c(rep('red', length(S10)),
                   rep('blue', length(S20))), pch=19,
             cex=1, cex.main=2, cex.axis=2,
             main = 't-SNE plot for RNA-seq and ATAC-seq',
             xlab='x', ylab='y')
        legend("topright", inset=.05,
               c("RNA-seq", "ATAC-seq"), fill=c('red', 'blue'), horiz=F, cex=2,
               bty = "n")
        dev.off()
    }
}
