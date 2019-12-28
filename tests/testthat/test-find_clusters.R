context("find_clusters")

test_that("ncores=2", {
  sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)

  sce1 <- scater::logNormCounts(sce1)
  metadata(sce1)$log.exprs.offset <- 1
  trend <- tech_trend(sce1, ncores=2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend, assay.type="logcounts", max.rank=100))
  expect_equal(ncol(reducedDim(sce1, "PCA")), 75)

<<<<<<< HEAD
  sce1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE)
  expect_setequal(as.vector(table(sce1$Cluster)), c(30, 32))
  expect_equal(levels(sce1$Cluster), c("clus_1", "clus_2"))
 
  # method="spinglass"
  sce1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE, method="spinglass")
  expect_setequal(as.vector(table(sce1$Cluster)), c(28, 22))
  expect_equal(levels(sce1$Cluster), c("clus_1", "clus_2"))	
=======
  expect_warning(sce1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE))
  expect_setequal(as.vector(table(sce1$Cluster)), c(32,30))
  
  # method="spinglass"
  expect_warning(sce1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE, method="spinglass"))
  expect_setequal(as.vector(table(sce1$Cluster)), c(28,22))	
>>>>>>> 60a08fcf9d9194eb35033649ae65cc86b1305e04
})

