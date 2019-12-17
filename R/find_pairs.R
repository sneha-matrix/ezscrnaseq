#' Cell cycle training step of the pair-based prediction
#'
#' Cell cycle pairs creation using \pkg{scran} \code{sandbag}.
#'
#' @param genes Ensembl gene IDs for genes in \code{sce}.
#' @param G1 vector of column number having G1 genes expression
#' @param S vector of column number having S genes expression
#' @param G2M vector of column number having G2M genes expression
#' @inheritParams qc_metrics
#' @inheritParams scran::sandbag
#' @return A list containing training gene pairs for G1, S, and G2M phase.
#' @export

find_pairs <- function(sce, G1, S, G2M, genes=NULL){
  if (is.null(genes)) genes <- rownames(sce)
  stopifnot(!is.null(G1), !is.null(S), !is.null(G2M), !is.null(rownames(sce)),
            length(which(rownames(sce) %in% genes)) > 0)

  sce2 <- sce[which(rownames(sce) %in% genes),, drop=FALSE]
  pairs <- scran::sandbag(sce2, list(G1=G1, S=S, G2M=G2M))
  return(pairs)
}

