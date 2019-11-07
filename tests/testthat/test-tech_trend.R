context("tech_trend")

test_that("default", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE))
  sce1 <- scater::normalize(sce1)
  trend <- tech_trend(sce1, span=1, plot=FALSE)
  expect_equal(round(trend(c(10, 20)), 5), c(0.00289, -0.06269))
})

test_that("ncores=2", {
  sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=2,plot=FALSE, verbose=FALSE))
  sce1 <- scater::normalize(sce1)
  trend <- tech_trend(sce1, ncores=2, plot=FALSE)
  expect_equal(round(trend(c(10, 20)), 5), c(0.00289, -0.06269))
})
