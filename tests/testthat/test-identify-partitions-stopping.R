test_that("split is rejected when one child fails minsize", {
  
  W <- rep(c(0, 0, 1, 1), 3)
  G <- rep(c(0, 1, 0, 1), 3)
  Y <- seq_len(12)
  
  X <- data.frame(X.1 = c(rep(0, 8), rep(1, 4)))
  
  p <- rep(0.5, 12)
  Ne_list <- replicate(12, integer(0), simplify = FALSE)
  Ne <- lengths(Ne_list)
  
  result <- identify_partitions_nct(
    method = "singular",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    depth = 2,
    minsize = 2,
    N = 12,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = rep(1, 4)
  )
  
  expect_equal(nrow(result), 1)
  expect_equal(result$TERMINAL[1], "LEAF")
})

test_that("split is rejected when missing covariates prevent a valid partition", {
  
  W <- rep(c(0, 0, 1, 1), 5)
  G <- rep(c(0, 1, 0, 1), 5)
  Y <- seq_len(20)
  
  X <- data.frame(X.1 = c(rep(0, 8), rep(1, 11), NA_real_))
  
  p <- rep(0.5, 20)
  Ne_list <- replicate(20, integer(0), simplify = FALSE)
  Ne <- lengths(Ne_list)
  
  result <- identify_partitions_nct(
    method = "singular",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    depth = 2,
    minsize = 1,
    N = 20,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = rep(1, 4)
  )
  
  expect_equal(nrow(result), 1)
  expect_equal(result$TERMINAL[1], "LEAF")
})

