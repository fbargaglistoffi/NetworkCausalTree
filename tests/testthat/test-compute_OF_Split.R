test_that("compute_OF_Split runs on tiny data and returns valid structure", {
  
  N <- 8
  W <- c(0,0,1,1, 0,0,1,1)
  G <- c(0,1,0,1, 0,1,0,1)
  Y <- c(1,2,3,4, 5,6,7,8)
  X <- matrix(c(1,1,1,1, 3,3,3,3), ncol=1)
  colnames(X) <- "X.1"
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5,8,7)
  p <- rep(0.5, N)
  
  pe <- compute_population_effects(N, W, G, Y, p, Ne)
  
  result <- compute_OF_Split(
    method = "singular",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = pe,
    total_variance = NULL, nleafs = 1
  )
  
  expect_length(result, 3)
  expect_true(is.numeric(as.numeric(result[1])))
  expect_true(!is.na(result[1]))
  expect_true(result[1] >= 0)
  expect_true(result[2] %in% X)
  expect_true(is.character(result[3]) && nchar(result[3]) > 0)
})