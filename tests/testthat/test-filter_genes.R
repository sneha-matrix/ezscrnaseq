context("filter_genes")

test_that("defult cutoff=0", {
  sce1 <- filter_genes(sce, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE)
  sce2 <- filter_genes(sce, cutoff=0, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)
})

test_that("cutoff=0.1", {
  sce1 <- filter_genes(sce, cutoff=0.1, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE)
  sce2 <- filter_genes(sce, cutoff=0.1, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)
})

test_that("defult ncore=1", {
  sce1 <- filter_genes(sce, , cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce2 <- filter_genes(sce, cutoff=0, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)
})

test_that("negative tests", {
  # negative cutoff	
  expect_error(sce1 <- filter_genes(sce, cutoff=-10.1, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE))
  # ncore = 0
  expect_error(sce1 <- filter_genes(sce, cutoff=-10.1, ncores=0, plot=FALSE, write=FALSE, verbose=FALSE))
  # logical tests for argument 
  expect_error(sce1 <- filter_genes(sce, cutoff=-10.1, ncores=1, plot=1))
  expect_error(sce1 <- filter_genes(sce, cutoff=-10.1, ncores=1, write=1))
  expect_error(sce1 <- filter_genes(sce, cutoff=-10.1, ncores=1, verbose=1))
})
test_that("truth table", {

  expect_warning(data("sc_example_counts"))
  sc_example_counts_shift <- sc_example_counts
  sc <- sc_example_counts_shift[1:200,1:10]
  sc[sc > 0]<-0
  sc[1:100, 1:5] <- 3
  sc[101:200, 6:10] <- 7
  sc[1:3,]<-0
  scVar <- sc
  scVar1 <- SingleCellExperiment(assays = list(counts = scVar))
  sce1 <- filter_genes(scVar1, cutoff=0, ncores=1, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_equal(dim(sce1)[1], 197)

})
