test_that("identify_partitions_nct runs without error and returns a tree dataframe", {
  set.seed(100)
  
  N <- 6
  W <- c(0,1,0,1,0,1)
  G <- c(0,0,1,1,0,1)
  Y <- 1 + W + G + rnorm(N, 0, 0.1)
  p <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)
  
  X <- matrix(1:6, ncol = 1)
  colnames(X) <- "X.1"
  
  population_effects <- rep(1,4)
  # for composite, they should add to 1, but for singular, normalization doesn't
  # matter
  alpha <- beta <- gamma <- delta <- 1
  depth <- 1
  minsize <- 1
  
  result <- identify_partitions_nct(
    method = "singular",
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    depth = depth, minsize = minsize,
    N = N, W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = population_effects
  )
  
  expect_s3_class(result, "data.frame")
  
  # >= 1 row
  expect_gte(nrow(result), 1)
  
  expect_true(all(c("NODE","OF","NOBS","FILTER","TERMINAL") %in% colnames(result)))
  
  expect_true(all(result$TERMINAL %in% c("SPLIT", "LEAF")))
  
  expect_true(all(is.numeric(result$OF) | is.na(result$OF)))
})