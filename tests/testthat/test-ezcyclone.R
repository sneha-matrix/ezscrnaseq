context("ezcyclone")

test_that("ncore=1", {
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("S","S", "G1"))

  mypairs <- find_pairs(sce, genes.list=genes.list, G1=G1, S=S, G2M=G2M)
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=mypairs)
  expect_equal(assignments$phases[1:3], c("S","G2M", "G1"))

})

test_that("ncore=2", {
  assignments <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("S","S", "G1"))

  mypairs <- find_pairs(sce, genes.list=genes.list, G1=G1, S=S, G2M=G2M)
  assignments <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=mypairs)
  expect_equal(assignments$phases[1:3], c("S","G2M", "G1"))
})
