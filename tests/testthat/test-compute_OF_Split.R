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

test_that("compute_OF_Split uses sample-size-weighted objectives", {
  
  n <- 12
  X <- matrix(c(rep(0, 8), rep(1, 4)), ncol = 1)
  colnames(X) <- "X.1"
  
  W <- c(0, 0, 1, 1, 0, 0, 1, 1,
         0, 1, 0, 1)
  G <- c(0, 1, 0, 1, 0, 1, 0, 1,
         0, 0, 1, 1)
  Y <- seq_len(n)
  
  p <- rep(0.5, n)
  Ne_list <- replicate(n, integer(0), simplify = FALSE)
  Ne <- lengths(Ne_list)
  
  result <- compute_OF_Split(
    method = "singular",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = rep(1, 4),
    total_variance = 1, nleafs = 2
  )
  
  left_of <- compute_OF_Value(
    "singular", 1, 0, 0, 0,
    N = 8, W = W[1:8], G = G[1:8], Y = Y[1:8],
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = 1:8,
    population_effects = rep(1, 4),
    total_variance = 1, nleafs = 2
  )
  
  right_of <- compute_OF_Value(
    "singular", 1, 0, 0, 0,
    N = 4, W = W[9:12], G = G[9:12], Y = Y[9:12],
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = 9:12,
    population_effects = rep(1, 4),
    total_variance = 1, nleafs = 2
  )
  
  expected <- (8 / 12) * left_of + (4 / 12) * right_of
  
  expect_equal(as.numeric(result["of"]), expected)
})
