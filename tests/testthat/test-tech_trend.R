context("tech_trend")

test_that("default", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE)
  sce2 <- scater::logNormCounts(sce1)
  metadata(sce2)$log.exprs.offset <- 1
  trend <- tech_trend(sce2, span=1, plot=FALSE)
  expect_equal(round(trend(c(10, 20)), 5), c(0.0028, -0.0963))
})

test_that("ncores=2", {
  sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=2,plot=FALSE, verbose=FALSE)
  sce2 <- scater::logNormCounts(sce1)
  metadata(sce2)$log.exprs.offset <- 1
  trend <- tech_trend(sce2, ncores=2, plot=FALSE)
  expect_equal(round(trend(c(10, 20)), 5), c(0.0028, -0.0963))
})







  
