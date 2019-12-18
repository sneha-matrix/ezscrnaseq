context("find_pairs")

test_that("ncore=1", {
  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M)
  genes2 <- row.names(sce)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes2)
  expect_equal(mypairs, mypairs2)

  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes)
  genes2 <- rbind(genes,genes)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes)
  expect_equal(mypairs, mypairs2)
})
