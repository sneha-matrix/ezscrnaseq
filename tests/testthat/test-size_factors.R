context("size_factors")

test_that("min.size=5", {
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=1,plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)
})
