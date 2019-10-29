#' Make a technical trend & decompose the gene-level variance
#'
#' Make a technical trend using \pkg{scran} \code{makeTechTrend} and decompose the gene-level variance using \pkg{scran} \code{decomposeVar}.
#'
#' @param assay_type A string specifying which assay values to use, e.g., "counts" or "logcounts".
#' @inheritParams qc_metrics
#' @inheritParams scran::makeTechTrend
#' @inheritParams scran::trendVar
#' @inheritParams scran::decomposeVar
#' @inheritParams stats::loess
#' @return A function accepting a mean log-expression as input and returning the variance of the log-expression as the output
#' @export

tech_trend <- function(sce, dispersion=0, span=0.4, block=NA, design=NA, assay_type="logcounts", ncores=1, prefix=NULL, plot=TRUE){

  # no spike ref: https://github.com/MarioniLab/scran/issues/7
  # according to scran::multiBlockVar(), use logcounts for tech trend
  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- SnowParam(workers=ncores, type=cl_type)
  register(bpstart(bp))
  var_fit_trend <- makeTechTrend(dispersion=dispersion, x=sce, BPPARAM=bp)

  var_fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=span), use.spikes=FALSE, assay.type=assay_type)
  var_out <- decomposeVar(sce, fit=var_fit, block=block, design=design, assay.type=assay_type, BPPARAM=bp)
  bpstop(bp)

  if (plot){
    grDevices::pdf(paste(c(prefix, "mean_variance_trend.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    graphics::plot(var_out$mean, var_out$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
    graphics::curve(var_fit_trend(var_out$mean), col="dodgerblue", add=TRUE, lwd=2)
    graphics::curve(var_fit$trend(var_out$mean), col="red", add=TRUE, lwd=2)
    graphics::legend("topright", legend=c("Technical noise", "All variance"), lty=1, lwd=2, col=c("dodgerblue", "red"))
  }

  return(var_fit_trend)
}
