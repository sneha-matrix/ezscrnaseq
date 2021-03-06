#' Make a technical trend & decompose the gene-level variance
#'
#' Make a technical trend using \pkg{scran} \code{makeTechTrend} and decompose the gene-level variance using
#' \pkg{scran} \code{modelGeneVar}.
#'
#' @param assay_type A string specifying which assay values to use, e.g., "counts" or "logcounts".
#' @inheritParams qc_metrics
#' @inheritParams scran::makeTechTrend
#' @inheritParams scran::modelGeneVar
#' @return A function accepting a mean log-expression as input and returning the variance of the log-expression as the output
#' @export

tech_trend <- function(sce, dispersion=0, assay_type="logcounts", block=NULL, design=NULL, ncores=1, prefix=NULL, plot=TRUE){

  stopifnot(ncores > 0, dispersion >=0, is.logical(plot))

  # no spike ref: https://github.com/MarioniLab/scran/issues/7
  # according to scran::multiBlockVar(), use logcounts for tech trend
  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(BiocParallel::bpstart(bp))
  var_fit_trend <- scran::makeTechTrend(dispersion=dispersion, x=sce, BPPARAM=bp)

  var_tot <- scran::modelGeneVar(sce, assay.type=assay_type, block=block, design=design)
  var_fit <- scran::fitTrendVar(means=var_tot$mean, vars=var_tot$total)

  #deprecated function
  #var_fit <- trendVar(sce, parametric=FALSE, loess.args=list(span=span), use.spikes=FALSE, assay.type=assay_type)
  #var_out <- decomposeVar(sce, fit=var_fit, block=block, design=design, assay.type=assay_type, BPPARAM=bp)

  BiocParallel::bpstop(bp)

  if (plot){
    grDevices::pdf(paste(c(prefix, "mean_variance_trend.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    graphics::plot(var_tot$mean, var_tot$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
    graphics::curve(var_fit_trend, col="dodgerblue", add=TRUE, lwd=2)
    graphics::curve(var_fit$trend(x), col="red", add=TRUE, lwd=2)
    graphics::legend("topright", legend=c("Technical noise", "All variance"), lty=1, lwd=2, col=c("dodgerblue", "red"))
  }
  return(var_fit_trend)
}
