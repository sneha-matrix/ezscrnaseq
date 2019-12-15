context("find_markers")

test_that("find_markers", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE))

  expect_warning(sce1 <- scater::normalize(sce1))
  expect_warning(trend <- tech_trend(sce1, plot=FALSE))
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend, assay.type="logcounts", max.rank=100))

  expect_warning(sce1 <- find_clusters(sce1, snn_k=5, plot=FALSE, verbose=FALSE))

  expect_warning(clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)),  
                                  fdr_cutoff=Inf, write=FALSE))
  expect_equal(round(clus_df["ENSG00000126249","ave.count"],3), 59.442)
  expect_equal(round(clus_df["ENSG00000222808","ave.count"],3), 16.249)

  expect_warning(clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)), 
                               fdr_cutoff=Inf, pval.type="all", write=FALSE))
  expect_equal(round(clus_df["ENSG00000126249","ave.count"],3), 59.442)
  expect_equal(round(clus_df["ENSG00000222808","ave.count"],3), 16.249)
})
