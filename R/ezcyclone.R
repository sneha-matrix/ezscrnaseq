#' Cell cycle phase classification
#'
#' Cell cycle phase classification using \pkg{scran} \code{cyclone}
#'
#' @param organism hsa or mmu for human or mouse genes.
#' @param gene.names Ensembl gene IDs for genes in \code{sce}.
#' @param ncores Number of cores.
#' @param seed Random seed.
#' @inheritParams qc_metrics
#' @inheritParams scran::cyclone
#' @return A list containing phases, scores, and normalized.scores .
#' @export

ezcyclone <- function(sce, organism="hsa", gene.names=rownames(sce), ncores=1, seed=100, iter=1000, min.iter=100, min.pairs=50, verbose=TRUE){

  if (organism=="hsa"){
    pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
  } else if(organism=="mmu"){
    pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  }

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- SnowParam(workers=ncores, type=cl_type)
  register(bpstart(bp))
  set.seed(seed)
  assignments <- cyclone(sce, pairs=pairs, gene.names=gene.names, iter=iter, min.iter=min.iter, min.pairs=min.pairs,  BPPARAM=bp, verbose=verbose)
  bpstop(bp)

  return(assignments)
}
