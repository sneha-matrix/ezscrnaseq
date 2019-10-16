#' Find cell clusters
#'
#' Find cell clusterd using \pkg{igraph}.
#'
#' @param seed Random seed.
#' @param snn_k The number of nearest neighbors to consider during graph construction.
#' @param "walktrap" or "spinglass" for finding communities in graphs via short random walks or a spin-glass model and simulated annealing.
#' @param min_member Minimal number of cluster members.
#' @inheritParams qc_metrics
#' @inheritParams scran::buildSNNGraph
#' @inheritParams igraph::cluster_walktrap
#' @inheritParams igraph::cluster_spinglass
#' @return A SingleCellExperiment object with cell cluster information.
#' @export

find_clusters <- function(sce, use_dimred="PCA", seed=100, snn_k=10, ncores=2, method="walktrap", steps=4, spins=25, min_member=20, prefix=NULL,
                          plot=TRUE, verbose=TRUE){

  set.seed(seed)

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- SnowParam(workers=ncores, type=cl_type)
  register(bpstart(bp))
  # snn
  snn_gr <- buildSNNGraph(sce, use.dimred=use_dimred, k=snn_k, BPPARAM=bp)
  bpstop(bp)

  # cluster
  if(method=="walktrap"){
    cluster_out <- igraph::cluster_walktrap(snn_gr, steps=steps)

  } else if(method=="spinglass"){
    cluster_out <- igraph::cluster_spinglass(snn_gr, spins=spins)

  }

  sce$Cluster <- factor(cluster_out$membership)
  sce$Cluster <- factor(paste0("clus_", sce$Cluster), levels=paste0("clus_", levels(sce$Cluster)))
  nc <- table(sce$Cluster)

  if (verbose){
    cat("Clusters found:\n")
    print(kable(nc))
  }


  # modularity score
  ms <- igraph::modularity(cluster_out)
  if (verbose) cat("\nModularity score: ", ms, "\n")

  # total weight between nodes
  mod_out <- clusterModularity(snn_gr, sce$Cluster, get.values=TRUE)
  ratio <- log2(mod_out$observed / mod_out$expected + 1)

  if (plot){
    pdf(paste(c(prefix, "clusters_total_weights.pdf"), collapse="_"))
    on.exit(dev.off())
    pheatmap(ratio, scale="none", cluster_rows=FALSE, cluster_cols=FALSE, color=colorRampPalette(c("white", "blue"))(100))
  }

  # rm small clusters
  num_small <- sum(nc < min_member)
  if(num_small > 0){
    if (verbose) cat("Removeing", num_small, "clusters that have cells less than", min_member, "\n")
    sce <- sce[, sce$Cluster %in% names(nc)[nc >=min_member]]
    sce$Cluster <- droplevels(sce$Cluster)
  }

  return(sce)
}
