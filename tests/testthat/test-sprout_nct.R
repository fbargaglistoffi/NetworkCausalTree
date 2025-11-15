test_that("sprout_nct builds tree on tiny sample and returns valid structure", {
  set.seed(100)
  
  N <- 6
  W <- c(0,1,0,1,0,1)
  G <- c(0,0,1,1,0,1)
  Y <- 1 + W + G + rnorm(N, 0, 0.05)
  p <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)
  
  X <- matrix(1:6, ncol = 1)
  colnames(X) <- "1"
  
  K <- c(1,1,2,2,3,3)
  
  sampled_clusters <- c(1,2)
  
  population_effects <- rep(1,4)
  alpha <- beta <- gamma <- delta <- 1
  depth <- 1
  minsize <- 1
  
  result <- sprout_nct(
    method = "singular",
    sampled_clusters = sampled_clusters,
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    depth = depth, minsize = minsize,
    N = N, W = W, G = G, Y = Y, X = X, K = K, p = p, Ne = Ne,
    population_effects = population_effects, Ne_list = Ne_list
  )
  
  
  expect_s3_class(result, "data.frame")
  
  expect_true(all(c("NODE","OF","NOBS","FILTER","TERMINAL") %in% colnames(result)))
  
  expect_gte(nrow(result), 1)
  
  expect_true(all(result$TERMINAL %in% c("SPLIT", "LEAF")))
})