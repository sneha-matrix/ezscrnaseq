#' Cell cycle training step of the pair-based prediction
#'
#' Cell cycle pairs creation using \pkg{scran} \code{sandbag}.
#'
#' @param genes Ensembl gene IDs for genes in \code{sce}.
#' @param G1 Vector of column indices of G1 genes.
#' @param S Vector of column indices of S genes.
#' @param G2M Vector of column indices of G2M genes.
#' @inheritParams qc_metrics
#' @inheritParams scran::sandbag
#' @details Genes in all three phases are required.
#' @return A list containing training gene pairs for G1, S, and G2M phase.
#' @export

find_pairs <- function(sce, G1, S, G2M, genes=NULL){
  if (is.null(genes)) genes <- rownames(sce)
  stopifnot(is.numeric(G1), is.numeric(G2M), is.numeric(S), !is.null(rownames(sce)),
            length(which(rownames(sce) %in% genes)) > 0)

  sce2 <- sce[which(rownames(sce) %in% genes),, drop=FALSE]
  pairs <- scran::sandbag(sce2, list(G1=G1, S=S, G2M=G2M))
  return(pairs)
}

