context("qc_metrics")

test_that("by_nmads", {
  expect_warning(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(2000, 70))
})

test_that("ncores=2", {
  expect_warning(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=2, plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(2000, 70))
})

test_that("not_by_nmads", {
  expect_warning(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(4e5,600,20), plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(2000, 26))
})

test_that("large_nmads", {
  expect_warning(expect_error(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(10,3,3), plot=FALSE, write=FALSE, verbose=FALSE)))
})

test_that("small counts", {
  expect_warning(expect_error(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(9,600,20), plot=FALSE, write=FALSE, verbose=FALSE)))
})
