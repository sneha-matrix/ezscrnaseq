context("find_pairs")

test_that("tests", {
  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M)
  genes2 <- row.names(sce)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes2)
  expect_equal(mypairs, mypairs2)

  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes)
  genes2 <- rbind(genes,genes)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes2)
  expect_equal(mypairs, mypairs2)

 #test for nrow(sce) <1000
 sce2 <- head(sce,999)
 expect_error(find_pairs(sce2, G1=G1, S=S, G2M=G2M))

 #test for empty pairs
 sce3 <- SingleCellExperiment::rbind(rep(head(sce),times=200))
 expect_error(find_pairs(sce3, G1=G1, S=S, G2M=G2M))
})
