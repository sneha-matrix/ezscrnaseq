context("find_pairs")

test_that("ncore=1", {
  mypairs <- find_pairs(sce, genes.list=genes.list, G1=G1, S=S, G2M=G2M)
  expect_equal(mypairs$G1[1,1], "LRG_396")
  expect_equal(mypairs$S[1,1], "ENSG00000262473")
  expect_equal(mypairs$G2M[1,1], "ENSG00000214511")
})  	
