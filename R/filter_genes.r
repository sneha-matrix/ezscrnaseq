#' Filter genes base on avarage counts
#'
#' Filter genes that have avarage counts fewer than \code{cutoff)}
#'
#' @param cutoff The cutoff of avarage counts for filtering.
#' @inheritParams qc_metrics
#' @inheritParams scater::calcAverage
#' @return A SingleCellExperiment object.
#' @export


filter_genes <- function(sce, cutoff=0, use_size_factors=FALSE, prefix=NULL, plot=TRUE, write=TRUE){

  rowData(sce)$ave.count <- calcAverage(sce, use_size_factors=use_size_factors)
  keep <- rowData(sce)$ave.count > cutoff
  n_keep <- sum(keep)
  cat("\nNumber of genes kept:", n_keep, "\n")

  if(plot){
    pdf(paste(c(prefix, "average_counts.pdf"), collapse="_"))
    on.exit(dev.off())

    hist(log10(rowData(sce)$ave.count), breaks=100, main="", col="grey", xlab=expression(Log[10]~"average count"))
    if(cutoff > min(rowData(sce)$ave.count)) abline(v=log10(cutoff), lty=2)
  }

  g2k <- as.data.frame(table(keep))
  if(write) write.csv(g2k, paste(c(prefix, "genes_to_keep.csv"), collapse="_"), row.names=FALSE)

  sce <- sce[keep, ]
  return(sce)
}
