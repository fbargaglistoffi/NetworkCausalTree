test_that("compute_OF_Value works for singular case on tiny example", {

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
  
  N <- length(Y)

  # not actually used but need for argument
  pe <- compute_population_effects(N, W, G, Y, p, Ne)

  alpha <- 1
  beta <- 0
  gamma <- 0
  delta <- 0
  
  result <- compute_OF_Value(
    method = "singular",
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    N = N, W = W, G = G, Y = Y, p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = pe,
    total_variance = NULL, nleafs = NULL
  )
  
  expected <- (EffTau1000(N, W, G, Y, p, Ne))^2
  
  expect_equal(result, expected, tolerance = 1e-8)
})

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
  # for composite, they should add to 1, but for singluar, normalization doesnt
  # matter
  alpha <- beta <- gamma <- delta <- 1
  depth <- 1
  minsize <- 1

  res <- identify_partitions_nct(
      method = "singular",
      alpha = alpha, beta = beta, gamma = gamma, delta = delta,
      depth = depth, minsize = minsize,
      N = N, W = W, G = G, Y = Y, X = X,
      p = p, Ne = Ne, Ne_list = Ne_list,
      population_effects = population_effects
    )

  expect_s3_class(res, "data.frame")

  # >= 1 row
  expect_gte(nrow(res), 1)

  expect_true(all(c("NODE","OF","NOBS","FILTER","TERMINAL") %in% colnames(res)))

  expect_true(all(res$TERMINAL %in% c("SPLIT", "LEAF")))

  expect_true(all(is.numeric(res$OF) | is.na(res$OF)))
})

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

  res <- sprout_nct(
      method = "singular",
      sampled_clusters = sampled_clusters,
      alpha = alpha, beta = beta, gamma = gamma, delta = delta,
      depth = depth, minsize = minsize,
      N = N, W = W, G = G, Y = Y, X = X, K = K, p = p, Ne = Ne,
      population_effects = population_effects, Ne_list = Ne_list
    )


  expect_s3_class(res, "data.frame")

  expect_true(all(c("NODE","OF","NOBS","FILTER","TERMINAL") %in% colnames(res)))

  expect_gte(nrow(res), 1)

  expect_true(all(res$TERMINAL %in% c("SPLIT", "LEAF")))
})

test_that("compute_effects_nct runs on tiny tree and returns valid output", {
  set.seed(100)

  N <- 6
  W <- c(0,1,0,1,0,1)
  G <- c(0,0,1,1,0,1)
  Y <- 1 + W + G + rnorm(N, 0, 0.05)
  p <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)
  
  X <- matrix(1:6, ncol = 1)
  colnames(X) <- "X.1"

  nct_partition <- data.frame(
    NODE = c(1,2),
    OF = c(0,0.2),
    NOBS = c(N, 3),
    FILTER = c(NA, "data_tree$X.1 >= 3"),
    TERMINAL = c("PARENT", "LEAF")
  )
  
  output <- "estimation"
  minsize <- 1
  
  res <- compute_effects_nct(
      output = output,
      nct_partition = nct_partition,
      N = N,
      W = W, G = G, Y = Y, X = X,
      Ne = Ne, Ne_list = Ne_list,
      p = p, minsize = minsize
    )

  expect_s3_class(res, "data.frame")

  expect_gte(nrow(res), 1)

  required_cols <- c("NODE","FILTER","TERMINAL","NOBS_TR","NOBS_EST",
                     "EFF1000_EST","SE1000_EST",
                     "EFF1101_EST","SE1101_EST",
                     "EFF1110_EST","SE1110_EST",
                     "EFF0100_EST","SE0100_EST")
  
  expect_true(all(required_cols %in% colnames(res)))

  expect_true(all(sapply(res$EFF1000_EST, is.numeric)))
  expect_true(all(sapply(res$SE1000_EST, is.numeric)))
})

