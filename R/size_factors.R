#' Calculate deconvolving size factors from cell pools
#'
#' Calculate deconvolving size factors from cell pools using \pkg{ scran} \code{computeSumFactors}
#'
#' @param group.col column name of the grouping variable on the samples in colData(sce).
#' @param seed Random seed.
#' @inheritParams qc_metrics
#' @inheritParams scran::quickCluster
#' @inheritParams scran::computeSumFactors
#' @return A SingleCellExperiment object with size factors.
#' @export

size_factors <- function(sce, min.size=200, max.size=3000, min.mean=0.1, group.col=NULL, method="igraph", seed=100, ncores=1, prefix=NULL, plot=TRUE, verbose=TRUE){

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- SnowParam(workers=ncores, type=cl_type)
  register(bpstart(bp))

  suppressWarnings(set.seed(seed = 100, sample.kind = "Rounding"))
  clusters <- quickCluster(sce, min.size=min.size, max.size=max.size, min.mean=min.mean, method=method, BPPARAM=bp)

  if (verbose){
    cat("\nNumber of cells in clusters:\n")
    print(kable(table(clusters)))
  }

  sce <- computeSumFactors(sce, cluster=clusters, max.cluster.size=max.size, min.mean=min.mean, positive=TRUE, BPPARAM=bp)
  bpstop(bp)

  if (verbose){
    cat("\nSummary of size factors:\n")
    print(summary(sizeFactors(sce)))
  }


  #rm zero size factors
  if (any(sizeFactors(sce)==0)){
    if (verbose) cat("\nRemoveing", sum(sizeFactors(sce)==0), "cells that have size factor of zero.\n")
    sce <- sce[, sizeFactors(sce) > 0]
  }

  # scatter
  if (plot){
    if (is.null(sce$total_counts)) stop("Total counts are not calculated yet.")

    grDevices::pdf(paste(c(prefix, "size_factor_scatter.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    if (!is.null(group.col)){
      graphics::plot(sizeFactors(sce), sce$total_counts/1e3, log="xy", ylab="Library size (thousands)", xlab="Size factor", main="Size factors from deconvolution",
           col=scales::alpha(as.numeric(as.factor(colData(sce)[, group.col]))+1, 0.3), pch=16)
      legd <- sort(unique(colData(sce)[, group.col]))
      graphics::legend("topleft", col=2:(length(legd)+1), pch=16, cex=1.2, legend=legd)
    } else {
      graphics::plot(sizeFactors(sce), sce$total_counts/1e3, log="xy", ylab="Library size (thousands)", xlab="Size factor", main="Size factors from deconvolution",
           col= scales::alpha(1, 0.3), pch=16)
    }
  }
  return(sce)
}
