#' Find marker genes for cell clusters
#'
#' Find marker genes for cell clusters using \pkg{scran} \code{findMarkers}
#'
#' @param annot gene annotation
#' @param fdr_cutoff FDR cutoff for top marker genes
#' @inheritParams qc_metrics
#' @inheritParams scran::findMarkers
#' @return A data.frame for the statistics and annotation of top marker geness
#' @export

find_markers <- function(sce, clusters, annot, block=NULL, design=NULL, pval.type="any", lfc=1, direction="up", assay_type="logcounts", ncores=1,
                         fdr_cutoff=0.25, prefix=NULL, write=TRUE){

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- SnowParam(workers=ncores, type=cl_type)
  register(bpstart(bp))
  markers <- findMarkers(sce, clusters=clusters, block=block, design=design, pval.type=pval.type, lfc=lfc, direction=direction, assay.type=assay_type, BPPARAM=bp)
  bpstop(bp)

  clus <- names(markers)
  if (pval.type=="any"){
    marker_sets <- lapply(seq_along(markers), function(i) tryCatch(data.frame(markers[[i]][markers[[i]][,3] < fdr_cutoff, 1:3], Cluster=clus[i]), error=function(e) NULL))
    suffix <- "simes"
  } else if(pval.type=="all"){
    marker_sets <- lapply(seq_along(markers), function(i) tryCatch(data.frame(markers[[i]][markers[[i]][,2] < fdr_cutoff, 1:2], Cluster=clus[i]), error=function(e) NULL))
    suffix <- "iut"
  }

  if (all(sapply(marker_sets, is.null))) stop (paste0("No markers have FDR < ", fdr_cutoff))

  marker_sets <- lapply(marker_sets, function(m) {m$ID <- rownames(m); m})
  marker_sets <- Reduce(rbind, marker_sets)
  marker_sets <- ezlimma::df_signif(marker_sets, 3)
  marker_sets <- data.frame(marker_sets, annot[marker_sets$ID, ])

  if (write) utils::write.csv(marker_sets, paste0(paste(c(prefix, "clusters_top_markers", suffix), collapse="_"), ".csv"), na="", row.names=FALSE)

  return(marker_sets)
}
