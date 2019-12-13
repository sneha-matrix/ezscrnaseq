context("find_clusters")

test_that("ncores=2", {
  expect_warning(sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE))
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE))

  expect_warning(sce1 <- scater::normalize(sce1))
  #trend <- tech_trend(sce1, ncores=2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  #expect_warning(sce1 <- denoisePCA(sce1, technical=trend, assay.type="logcounts", max.rank=100))
  #expect_equal(ncol(reducedDim(sce1, "PCA")), 75)

  #expect_warning(sce1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE))
  #expect_equal(as.vector(table(sce1$Cluster)), c(30,32))
})

