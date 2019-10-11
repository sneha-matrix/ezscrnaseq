context("qc_metrics")

test_that("by_nmads", {
  sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), plot=FALSE, write=FALSE)
  expect_equal(dim(sce1), c(2000, 34))
})

test_that("not_by_nmads", {
  expect_warning(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(4e5,600,20), plot=FALSE, write=FALSE))
  expect_equal(dim(sce1), c(2000, 13))
})
