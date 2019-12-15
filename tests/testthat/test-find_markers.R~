context("find_markers")

test_that("find_markers", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  expect_warning(sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE))

  sce1 <- scater::normalize(sce1)
  trend <- tech_trend(sce1, plot=FALSE)
  set.seed(seed = 100, sample.kind = "Rounding")
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend, assay.type="logcounts", max.rank=100))

  sce1 <- find_clusters(sce1, snn_k=5, plot=FALSE, verbose=FALSE)

  clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("ENSG00000227430", "ENSG00000277058", "ENSG00000241846"))

  clus_df <- find_markers(sce1, clusters=sce1$Cluster, annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, pval.type="all", write=FALSE)
  expect_equal(rownames(clus_df)[1:3], c("ENSG00000227430", "ENSG00000277058", "ENSG00000241846"))
})
