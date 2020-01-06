context("find_markers")

test_that("find_markers", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE)

  sce2 <- scater::logNormCounts(sce1)
  metadata(sce2)$log.exprs.offset <- 1
  trend <- tech_trend(sce2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce3 <- denoisePCA(sce2, technical=trend, assay.type="logcounts", max.rank=100))

  sceCluster3 <- find_clusters(sce3, snn_k=5, plot=FALSE, verbose=FALSE)

  clus_df <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1), 
				row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)
  expect_equal(levels(factor(clus_df$Cluster)), c("clus_1", "clus_2"))
  expect_equal(round(clus_df["ENSG00000126249","ave.count"],3), 59.442)
  expect_equal(round(clus_df["ENSG00000222808","ave.count"],3), 16.249)
})


