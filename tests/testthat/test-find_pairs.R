context("find_pairs")

test_that("ncore=1", {
  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes.list=genes.list)
  expect_equal(mypairs$G1[1,1], "LRG_396")
  expect_equal(mypairs$S[1,1], "ENSG00000262473")
  expect_equal(mypairs$G2M[1,1], "ENSG00000214511")

  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M)
  expect_equal(mypairs$G1[1,1], "ENSG00000271537")
  expect_equal(mypairs$S[1,1], "ENSG00000214511")
  expect_equal(mypairs$G2M[1,1], "ENSG00000214511")
  
  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M)
  genes.list2 <- row.names(sce)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes.list=genes.list2)
  expect_equal(mypairs, mypairs2)

  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes.list=genes.list)
  genes.list2 <- rbind(genes.list,genes.list)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes.list=genes.list)
  expect_equal(mypairs, mypairs2)
})  	
