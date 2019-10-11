context("ezcyclone")

test_that("ncore=1", {
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("G1","G1","G2M"))
  #expect_equal(assignments$scores[1, 2], 0.615)
})

test_that("ncore=2", {
  assignments <- ezcyclone(sce, organism="hsa", ncore=2, min.pairs=5, verbose=FALSE, iter=10000, min.iter=1000)
  expect_equal(assignments$phases[1:3], c("G1","G1","G2M"))
  #expect_equal(assignments$scores[1, 2], 0.597)
})
