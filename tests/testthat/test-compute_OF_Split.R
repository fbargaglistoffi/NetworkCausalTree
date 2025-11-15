test_that("compute_OF_Split runs on tiny data and returns valid structure", {
  
  Ne_list <- list(
    c(2),
    c(1,3),
    c(2)
  )
  
  Ne <- c(1,2,1)
  p  <- c(0.5, 0.8, 0.2)
  
  W <- c(1,0,1)
  G <- c(0,1,0)
  
  Y <- c(10, 5, 15)
  
  X <- matrix(c(1,2,3), ncol=1)
  colnames(X) <- "X.1"
  
  N <- length(Y)
  
  pe <- compute_population_effects(N, W, G, Y, p, Ne)
  
  result <- compute_OF_Split(
    method = "singular",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    N = N, W = W, G = G, Y = Y, X = X,
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