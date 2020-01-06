context("size_factors")

test_that("default vs user provided args", {
  #test for ncore	
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for min.size = 10	
  sce1 <- size_factors(sce, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=10, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   max.size = 3000	
  sce1 <- size_factors(sce, min.size=5, min.mean=0, max.size = 3000, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   min.mean = 0.1	
  sce1 <- size_factors(sce, min.size=5, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0.1 , ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   seed = 100	
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, seed = 100, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   method = igraph 
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, method="igraph", plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #negative test for min.size > max.size	
  expect_error(size_factors(sce, min.size=500, max.size=5, min.mean=0, ncores=1, method="mycluster", plot=FALSE))
  # min.mean < 0
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=1 , seed=1 , min.mean= -1, plot=FALSE))
  # ncores = 0
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=0 , seed=1 , min.mean= 0, plot=FALSE))
  # seed =0
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=1 , seed=0 , min.mean= 1, plot=FALSE))

})



