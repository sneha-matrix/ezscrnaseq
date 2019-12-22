context("ezcyclone")

test_that("ncore=1", {
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("S","S", "G1"))

  mypairs <- find_pairs(sce, genes=genes, G1=G1, S=S, G2M=G2M)
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=mypairs)
  expect_equal(assignments$phases[1:3], c("S","G2M", "G1"))

})

test_that("ncore=2", {
  assignments <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("S","S", "G1"))

  mypairs <- find_pairs(sce, genes=genes, G1=G1, S=S, G2M=G2M)
  assignments <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=mypairs)
  expect_equal(assignments$phases[1:3], c("S","G2M", "G1"))
})

test_that("negative tests",{
  emptypairs <- scran::sandbag(head(sce), list(G1=G1, S=S, G2M=G2M))
  expect_error(ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=emptypairs))
  
  #pairs which return NA
  p2 <- list()		
  p2$S <- as.data.frame(t(genes[3:4]))
  colnames(p2$S)<-c("first","second")
  p2$G1 <- as.data.frame(t(genes[6:7]))
  colnames(p2$G1)<-c("first","second")
  p2$G2M <- as.data.frame(t(genes[1:2]))
  colnames(p2$G2M)<-c("first","second")
  expect_error(ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000, pairs=p2))
})
