test_that("NetworkCausalTree runs and returns expected structure", {
  set.seed(100)
  
  N <- 80
  M <- 2
  k <- 4
  
  K <- rep(1:k, each = N / k)
  
  X <- matrix(
    rbinom(N * M, size = 1, prob = 0.5),
    nrow = N,
    ncol = M
  )
  
  A <- matrix(0, nrow = N, ncol = N)
  
  for (cluster in unique(K)) {
    idx <- which(K == cluster)
    
    for (j in seq_along(idx)) {
      current <- idx[j]
      next_unit <- idx[ifelse(j == length(idx), 1, j + 1)]
      
      A[current, next_unit] <- 1
      A[next_unit, current] <- 1
    }
  }
  
  treatment_pattern <- c(0, 0, 0, 1, 0, 1, 1, 0, 0, 0)
  W <- rep(treatment_pattern, length.out = N)
  
  Y <- rnorm(N)
  
  p <- rep(0.5, N)
  
  result <- NetworkCausalTree(
    X = X,
    Y = Y,
    W = W,
    A = A,
    K = K,
    p = p,
    effect_weights = c(1, 0, 0, 0),
    ratio_disc = 0.5,
    depth = 3,
    minsize = 2,
    method = "singular",
    output = "estimation"
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("FILTER", "NOBS_EST") %in% names(result)))
  expect_gt(nrow(result), 0)
})

test_that("NOBS_EST at root matches held-out estimation sample size", {
  
  set.seed(100)
  d <- data_generator_direct_indirect(
    N = 500,
    k = 10,
    M = 2,
    p = rep(0.3, 500),
    remove_isolates = TRUE
  )
  
  ratio_disc <- 0.5
  k <- length(unique(d$K))
  nclusters_disc <- round(k * ratio_disc)
  
  set.seed(200)
  clusters_disc <- sample(
    unique(d$K),
    size = nclusters_disc,
    replace = FALSE
  )
  
  clusters_est <- setdiff(unique(d$K), clusters_disc)
  expected_nobs_est <- sum(d$K %in% clusters_est)
  
  set.seed(200)
  fit <- NetworkCausalTree(
    Y = d$Y,
    W = d$W,
    X = d$X,
    A = d$A,
    K = d$K,
    p = d$p,
    effect_weights = c(1, 0, 0, 0),
    method = "singular",
    depth = 1,
    minsize = 20,
    ratio_disc = ratio_disc
  )
  
  root_idx <- which(fit$FILTER == "NA")
  
  expect_length(root_idx, 1)
  expect_equal(fit$NOBS_EST[root_idx], expected_nobs_est)
})

test_that("held-out estimation outcomes do not affect composite splitting", {
  set.seed(100)
  
  N <- 400
  k <- 10
  cluster_size <- N / k
  
  K <- rep(1:k, each = cluster_size)
  
  X <- matrix(
    rep(rep(c(0, 1), each = 20), k),
    ncol = 1
  )
  
  A <- matrix(0, N, N)
  
  for (cluster in unique(K)) {
    idx <- which(K == cluster)
    
    for (j in seq_along(idx)) {
      current <- idx[j]
      next_unit <- idx[ifelse(j == length(idx), 1, j + 1)]
      
      A[current, next_unit] <- 1
      A[next_unit, current] <- 1
    }
  }
  
  W <- rep(c(0, 0, 0, 1, 1, 1, 0, 1), length.out = N)
  
  G <- as.numeric(as.vector(A %*% W) > 0)
  
  p <- rep(0.3, N)
  
  Y <- 2 +
    (1 + 5 * X[, 1]) * W +
    (1 + 3 * X[, 1]) * W * G +
    rnorm(N, 0, 0.01)
  
  ratio_disc <- 0.5
  
  set.seed(200)
  clusters_disc <- sample(unique(K), size = round(k * ratio_disc))
  
  clusters_est <- setdiff(unique(K), clusters_disc)
  est_units <- which(K %in% clusters_est)
  
  Y_modified <- Y
  Y_modified[est_units] <- Y_modified[est_units] + 10000
  
  set.seed(200)
  fit_original <- NetworkCausalTree(
    Y = Y, W = W, X = X, A = A, K = K, p = p,
    effect_weights = c(0.5, 0.5, 0, 0),
    method = "composite",
    depth = 1,
    minsize = 2,
    ratio_disc = ratio_disc,
    output = "detection"
  )
  
  set.seed(200)
  fit_modified <- NetworkCausalTree(
    Y = Y_modified, W = W, X = X, A = A, K = K, p = p,
    effect_weights = c(0.5, 0.5, 0, 0),
    method = "composite",
    depth = 1,
    minsize = 2,
    ratio_disc = ratio_disc,
    output = "detection"
  )
  
  expect_gt(nrow(fit_original), 1)
  
  expect_equal(fit_original$FILTER, fit_modified$FILTER)
  expect_equal(fit_original$NOBS_TR, fit_modified$NOBS_TR)
  expect_equal(fit_original$OF, fit_modified$OF)
})
