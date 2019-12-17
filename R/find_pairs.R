#' Cell cycle training step of the pair-based prediction  
#'
#' Cell cycle pairs creation using \pkg{scran} \code{sandbag}
#'
#' @param sce A numeric matrix-like object of gene expression values where rows are genes and columns are cells.
#' Alternatively, a SingleCellExperiment object containing such a matrix.
#' @param genes.list Ensembl gene IDs for genes in \code{sce}.
#' @param G1 vector of column number having G1 genes expression
#' @param S vector of column number having S genes expression
#' @param G2M vector of column number having G2M genes expression
#' @inheritParams scran::sandbag
#' @return A list containing training gene pairs for G1, S and G2M phase.
#' @export

find_pairs <- function(sce, G1, S, G2M, genes.list=NULL){
  stopifnot(!is.null(G1), !is.null(S), !is.null(G2M))
  if(is.null(genes.list)) genes.list=row.names(sce)

  sce2 <- sce[row.names(sce) %in% genes.list,]
  pairs <- sandbag(sce2, list(G1=G1, S=S, G2M=G2M))
  return(pairs)
}

