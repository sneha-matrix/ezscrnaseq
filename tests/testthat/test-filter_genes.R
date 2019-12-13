context("filter_genes")

test_that("cutoff=0", {
  expect_warning(sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(1980, 80))
})

test_that("ncores=2", {
  expect_warning(sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(1980, 80))
})

test_that("cutoff=0.1", {
  expect_warning(sce1 <- filter_genes(sce, cutoff=0.1, plot=FALSE, write=FALSE, verbose=FALSE))
  expect_equal(dim(sce1), c(1943, 80))
})
