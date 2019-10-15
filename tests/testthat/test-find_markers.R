context("find_markers")

test_that("clusters", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE))

  sce1 <- scater::normalize(sce1)
  trend <- tech_trend(sce1, plot=FALSE)
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend, approximate=TRUE, rand.seed=100, assay.type="logcounts", max.rank=100))
  expect_equal(ncol(reducedDim(sce1, "PCA")), 74)

  sce1 <- find_clusters(sce1, snn_k=5, plot=FALSE, verbose=FALSE)

  clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf)
  expect_equal(rownames(clus_df)[1:3], c("ENSG00000268181", "ENSG00000182583", "ENSG00000207065"))

  clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, pval.type="all")
  expect_equal(rownames(clus_df)[1:3], c("ENSG00000268181", "ENSG00000182583", "ENSG00000207065"))
})
