context("filter_genes")

test_that("cutoff=0", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE)
  expect_equal(dim(sce1), c(1955, 40))
})

test_that("cutoff=0.1", {
  sce1 <- filter_genes(sce, cutoff=0.1, plot=FALSE, write=FALSE)
  expect_equal(dim(sce1), c(1919, 40))
})
