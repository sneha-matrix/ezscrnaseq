context("size_factors")

test_that("min.size=5", {
  expect_warning(sce1 <- size_factors(sce, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE))
  expect_equal(round(sizeFactors(sce1)[1:3], 3), c(0.459, 1.109, 1.057))
})

test_that("ncores=2", {
  expect_warning(sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE))
  expect_equal(round(sizeFactors(sce1)[1:3], 3), c(0.459, 1.109, 1.057))
})
