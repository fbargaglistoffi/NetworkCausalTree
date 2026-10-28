test_that("pi preserves neighbors outside a tree leaf", {
  p <- c(0.2, 0.5, 0.8)
  
  Ne_list <- list(
    c(2, 3),
    c(1),
    c(1)
  )
  
  leaf_ids <- 1
  
  result <- pi(
    i = leaf_ids,
    w = 1,
    g = 0,
    p = p,
    Ne_list = Ne_list
  )
  
  expected <- 0.2 * (1 - 0.5) * (1 - 0.8)
  
  expect_equal(result, expected)
})

test_that("EffTau1000 preserves full-network probabilities in a leaf", {
  p <- c(0.2, 0.5, 0.8)
  
  Ne_list <- list(
    c(2, 3),
    c(1),
    c(1)
  )
  
  W <- c(1, 0)
  G <- c(0, 0)
  Y <- c(2, 3)
  
  result <- EffTau1000(
    N = 2, W = W, G = G, Y = Y,
    p = p, Ne_list = Ne_list,
    unit_ids = c(1, 2)
  )
  
  expected <- (2 / 0.02 - 3 / 0.4) / 2
  
  expect_equal(result, expected)
})

test_that("EffTau1000 preserves probabilities after subsetting", {
  p <- c(0.2, 0.5, 0.8)
  
  Ne_list <- list(
    c(2, 3),
    c(1),
    c(1)
  )
  
  W <- c(1, 0, 0)
  G <- c(0, 0, 0)
  Y <- c(2, 3, 0)
  
  full_result <- EffTau1000(
    N = 3, W = W, G = G, Y = Y,
    p = p, Ne_list = Ne_list
  )
  
  leaf_result <- EffTau1000(
    N = 2, W = W[1:2], G = G[1:2], Y = Y[1:2],
    p = p, Ne_list = Ne_list,
    unit_ids = c(1, 2)
  )
  
  expect_equal(full_result * 3, leaf_result * 2)
})

test_that("pij accounts for shared neighbors", {
  p <- rep(0.5, 3)
  
  Ne_list <- list(
    3,
    3,
    c(1, 2)
  )
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 0,
    wj = 1, gj = 0,
    p = p, Ne_list = Ne_list
  )
  
  expected <- 0.5 * 0.5 * 0.5
  
  expect_equal(result, expected)
})

test_that("pij returns zero for incompatible exposures", {
  p <- c(0.5, 0.5)
  
  Ne_list <- list(2, 1)
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 0,
    wj = 1, gj = 0,
    p = p, Ne_list = Ne_list
  )
  
  expect_equal(result, 0)
})

test_that("pij handles heterogeneous assignment probabilities", {
  p <- c(0.2, 0.6, 0.8)
  
  Ne_list <- list(
    3,
    3,
    c(1, 2)
  )
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 0,
    wj = 0, gj = 0,
    p = p, Ne_list = Ne_list
  )
  
  expected <- 0.2 * (1 - 0.6) * (1 - 0.8)
  
  expect_equal(result, expected)
})

test_that("pij agrees with exact enumeration", {
  p <- c(0.2, 0.6, 0.8)
  
  Ne_list <- list(
    c(2, 3),
    c(1, 3),
    c(1, 2)
  )
  
  assignments <- expand.grid(
    W1 = 0:1, W2 = 0:1, W3 = 0:1
  )
  
  for (wi in 0:1) {
    for (gi in 0:1) {
      for (wj in 0:1) {
        for (gj in 0:1) {
          
          exact <- 0
          
          for (k in seq_len(nrow(assignments))) {
            W <- as.numeric(assignments[k, ])
            
            G1 <- as.integer(any(W[Ne_list[[1]]] == 1))
            G2 <- as.integer(any(W[Ne_list[[2]]] == 1))
            
            if (W[1] == wi && G1 == gi &&
                W[2] == wj && G2 == gj) {
              
              exact <- exact + prod(
                ifelse(W == 1, p, 1 - p)
              )
            }
          }
          
          result <- pij(
            i = 1, j = 2,
            wi = wi, gi = gi,
            wj = wj, gj = gj,
            p = p, Ne_list = Ne_list
          )
          
          expect_equal(result, exact, tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("Vary handles zero joint exposure probabilities", {
  p <- rep(0.5, 3)
  
  Ne_list <- list(2, c(1, 3), 2)
  Ne <- c(1, 2, 1)
  
  W <- c(1, 0, 1)
  G <- c(0, 1, 0)
  Y <- c(2, 3, 4)
  
  result <- Vary(
    w = 1, g = 0,
    Y = Y, W = W, G = G,
    p = p, Ne = Ne, Ne_list = Ne_list
  )
  
  expected <- 240 + 128 + 80
  
  expect_equal(result, expected)
})

test_that("Covy applies zero-probability correction correctly", {
  p <- rep(0.5, 8)
  
  Ne_list <- list(
    2, 1,
    4, 3,
    6, 5,
    8, 7
  )
  
  Ne <- rep(1, 8)
  
  W <- c(1, 0, 1, 0, 0, 0, 0, 0)
  G <- c(0, 1, 0, 1, 0, 0, 0, 0)
  Y <- rep(1, 8)
  
  result <- Covy(
    N = 8,
    w1 = 1, g1 = 0,
    w2 = 0, g2 = 0,
    W = W, G = G, Y = Y,
    p = p, Ne = Ne, Ne_list = Ne_list
  )
  
  expect_equal(result, -24)
})

test_that("Covy is conservative under exact randomization", {
  
  p <- c(0.5, 0.5)
  Ne_list <- list(2, 1)
  Ne <- c(1, 1)
  
  Y00 <- c(2, 3)
  Y10 <- c(4, 5)
  
  assignments <- expand.grid(W1 = 0:1, W2 = 0:1)
  
  ht_a <- numeric(nrow(assignments))
  ht_b <- numeric(nrow(assignments))
  cov_est <- numeric(nrow(assignments))
  
  for (k in seq_len(nrow(assignments))) {
    
    W <- as.numeric(assignments[k, ])
    G <- c(W[2], W[1])
    
    Y <- numeric(2)
    
    for (i in 1:2) {
      if (W[i] == 0 && G[i] == 0) {
        Y[i] <- Y00[i]
      } else if (W[i] == 1 && G[i] == 0) {
        Y[i] <- Y10[i]
      } else {
        Y[i] <- 0
      }
    }
    
    ht_a[k] <- sum(
      (W == 0 & G == 0) * Y /
        pi(1:2, 0, 0, p, Ne_list)
    ) / 2
    
    ht_b[k] <- sum(
      (W == 1 & G == 0) * Y /
        pi(1:2, 1, 0, p, Ne_list)
    ) / 2
    
    cov_est[k] <- Covy(
      N = 2,
      w1 = 0, g1 = 0,
      w2 = 1, g2 = 0,
      W = W, G = G, Y = Y,
      p = p, Ne = Ne, Ne_list = Ne_list
    ) / 2^2
  }
  
  true_cov <- cov(ht_a, ht_b) * 3 / 4
  expected_cov_est <- mean(cov_est)
  
  expect_lte(expected_cov_est, true_cov)
})

test_that("Vary is conservative under exact randomization", {
  
  p <- c(0.5, 0.5)
  Ne_list <- list(2, 1)
  Ne <- c(1, 1)
  
  Y00 <- c(2, 3)
  
  assignments <- expand.grid(W1 = 0:1, W2 = 0:1)
  
  ht <- numeric(nrow(assignments))
  var_est <- numeric(nrow(assignments))
  
  for (k in seq_len(nrow(assignments))) {
    
    W <- as.numeric(assignments[k, ])
    G <- c(W[2], W[1])
    
    Y <- ifelse(W == 0 & G == 0, Y00, 0)
    
    ht[k] <- sum(
      (W == 0 & G == 0) * Y /
        pi(1:2, 0, 0, p, Ne_list)
    ) / 2
    
    var_est[k] <- Vary(
      w = 0, g = 0,
      Y = Y, W = W, G = G,
      p = p, Ne = Ne, Ne_list = Ne_list
    ) / 2^2
  }
  
  true_var <- mean((ht - mean(ht))^2)
  expected_var_est <- mean(var_est)
  
  expect_gte(expected_var_est, true_var)
})

test_that("Vartau1000 is conservative under exact randomization", {
  
  p <- c(0.5, 0.5)
  Ne_list <- list(2, 1)
  Ne <- c(1, 1)
  
  Y00 <- c(2, 3)
  Y10 <- c(4, 5)
  
  assignments <- expand.grid(W1 = 0:1, W2 = 0:1)
  
  tau_hat <- numeric(nrow(assignments))
  var_est <- numeric(nrow(assignments))
  
  for (k in seq_len(nrow(assignments))) {
    
    W <- as.numeric(assignments[k, ])
    G <- c(W[2], W[1])
    
    Y <- numeric(2)
    
    for (i in seq_along(W)) {
      if (W[i] == 0 && G[i] == 0) {
        Y[i] <- Y00[i]
      } else if (W[i] == 1 && G[i] == 0) {
        Y[i] <- Y10[i]
      } else {
        Y[i] <- 0
      }
    }
    
    pi00 <- pi(1:2, 0, 0, p, Ne_list)
    pi10 <- pi(1:2, 1, 0, p, Ne_list)
    
    ht00 <- sum((W == 0 & G == 0) * Y / pi00) / 2
    ht10 <- sum((W == 1 & G == 0) * Y / pi10) / 2
    
    tau_hat[k] <- ht10 - ht00
    
    var_est[k] <- Vartau1000(
      N = 2, Y = Y, W = W, G = G,
      p = p, Ne = Ne, Ne_list = Ne_list
    )
  }
  
  true_var <- mean((tau_hat - mean(tau_hat))^2)
  expected_var <- mean(var_est)
  
  expect_true(all(is.finite(var_est)))
  expect_gte(expected_var, true_var)
})

test_that("Vartau1000 is conservative on a three-unit network", {
  
  # Star network: units 2 and 3 share neighbor 1
  Ne_list <- list(c(2, 3), 1, 1)
  Ne <- c(2, 1, 1)
  p <- c(0.3, 0.5, 0.7)
  
  # Fixed potential outcomes for each exposure condition
  Y00 <- c(2, 3, 4)
  Y10 <- c(5, 6, 7)
  Y01 <- c(3, 4, 5)
  Y11 <- c(6, 7, 8)
  
  assignments <- expand.grid(W1 = 0:1, W2 = 0:1, W3 = 0:1)
  
  tau_hat <- numeric(nrow(assignments))
  var_est <- numeric(nrow(assignments))
  assignment_prob <- numeric(nrow(assignments))
  
  for (k in seq_len(nrow(assignments))) {
    
    W <- as.numeric(assignments[k, ])
    
    G <- c(
      as.integer(any(W[c(2, 3)] == 1)),
      W[1],
      W[1]
    )
    
    Y <- numeric(3)
    
    for (i in seq_along(W)) {
      if (W[i] == 0 && G[i] == 0) Y[i] <- Y00[i]
      if (W[i] == 1 && G[i] == 0) Y[i] <- Y10[i]
      if (W[i] == 0 && G[i] == 1) Y[i] <- Y01[i]
      if (W[i] == 1 && G[i] == 1) Y[i] <- Y11[i]
    }
    
    pi00 <- pi(1:3, 0, 0, p, Ne_list)
    pi10 <- pi(1:3, 1, 0, p, Ne_list)
    
    ht00 <- sum((W == 0 & G == 0) * Y / pi00) / 3
    ht10 <- sum((W == 1 & G == 0) * Y / pi10) / 3
    
    tau_hat[k] <- ht10 - ht00
    
    var_est[k] <- Vartau1000(
      N = 3, Y = Y, W = W, G = G,
      p = p, Ne = Ne, Ne_list = Ne_list
    )
    
    assignment_prob[k] <- prod(p^W * (1 - p)^(1 - W))
  }
  
  true_mean <- sum(assignment_prob * tau_hat)
  
  true_var <- sum(
    assignment_prob * (tau_hat - true_mean)^2
  )
  
  expected_var <- sum(assignment_prob * var_est)
  
  expect_true(all(is.finite(var_est)))
  expect_equal(sum(assignment_prob), 1)
  expect_gte(expected_var + 1e-10, true_var)
})

test_that("Covy is conservative on a three-unit network", {
  
  Ne_list <- list(c(2, 3), 1, 1)
  Ne <- c(2, 1, 1)
  p <- c(0.3, 0.5, 0.7)
  
  Y00 <- c(2, 3, 4)
  Y10 <- c(5, 6, 7)
  Y01 <- c(3, 4, 5)
  Y11 <- c(6, 7, 8)
  
  assignments <- expand.grid(
    W1 = 0:1, W2 = 0:1, W3 = 0:1
  )
  
  ht00 <- numeric(nrow(assignments))
  ht10 <- numeric(nrow(assignments))
  cov_est <- numeric(nrow(assignments))
  assignment_prob <- numeric(nrow(assignments))
  
  for (k in seq_len(nrow(assignments))) {
    
    W <- as.numeric(assignments[k, ])
    
    G <- c(
      as.integer(any(W[c(2, 3)] == 1)),
      W[1],
      W[1]
    )
    
    Y <- numeric(3)
    
    for (i in seq_along(W)) {
      if (W[i] == 0 && G[i] == 0) Y[i] <- Y00[i]
      if (W[i] == 1 && G[i] == 0) Y[i] <- Y10[i]
      if (W[i] == 0 && G[i] == 1) Y[i] <- Y01[i]
      if (W[i] == 1 && G[i] == 1) Y[i] <- Y11[i]
    }
    
    ht00[k] <- sum(
      (W == 0 & G == 0) * Y /
        pi(1:3, 0, 0, p, Ne_list)
    ) / 3
    
    ht10[k] <- sum(
      (W == 1 & G == 0) * Y /
        pi(1:3, 1, 0, p, Ne_list)
    ) / 3
    
    cov_est[k] <- Covy(
      N = 3,
      w1 = 0, g1 = 0,
      w2 = 1, g2 = 0,
      W = W, G = G, Y = Y,
      p = p, Ne = Ne, Ne_list = Ne_list
    ) / 3^2
    
    assignment_prob[k] <- prod(
      p^W * (1 - p)^(1 - W)
    )
  }
  
  mean00 <- sum(assignment_prob * ht00)
  mean10 <- sum(assignment_prob * ht10)
  
  true_cov <- sum(
    assignment_prob *
      (ht00 - mean00) *
      (ht10 - mean10)
  )
  
  expected_cov <- sum(assignment_prob * cov_est)
  
  expect_true(all(is.finite(cov_est)))
  expect_equal(sum(assignment_prob), 1)
  expect_lte(expected_cov, true_cov + 1e-10)
})

test_that("pij detects structural zeros for identical neighborhoods", {
  
  Ne_list <- list(3, 3, c(1, 2))
  p <- c(0.3, 0.6, 0.8)
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 0,
    wj = 0, gj = 1,
    Ne_list = Ne_list, p = p
  )
  
  expect_identical(result, 0)
})

test_that("pij preserves genuinely small positive probabilities", {
  
  Ne_list <- list(
    c(3, 4),
    3,
    c(1, 2),
    1
  )
  
  p <- c(0.5, 0.5, 0.8, 1e-10)
  
  result <- pij(
    i = 1, j = 2,
    wi = 0, gi = 1,
    wj = 0, gj = 0,
    Ne_list = Ne_list, p = p
  )
  
  expected <- 0.5 * 0.5 * (1 - 0.8) * 1e-10
  
  expect_gt(result, 0)
  expect_equal(result, expected, tolerance = 1e-8)
})

test_that("Vary and Covy preserve neighbors outside the leaf", {
  
  Ne_list <- list(
    3,
    3,
    c(1, 2)
  )
  
  p <- c(0.5, 0.5, 0.5)
  Ne <- lengths(Ne_list)
  
  W_full <- c(0, 1, 0)
  G_full <- c(0, 0, 1)
  Y_full <- c(2, 4, 6)
  
  unit_ids <- c(1, 2)
  
  W <- W_full[unit_ids]
  G <- G_full[unit_ids]
  Y <- Y_full[unit_ids]
  
  vary_result <- Vary(
    w = 0, g = 0,
    Y = Y, W = W, G = G,
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = unit_ids
  )
  
  covy_result <- Covy(
    N = length(unit_ids),
    w1 = 0, g1 = 0,
    w2 = 1, g2 = 0,
    W = W, G = G, Y = Y,
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = unit_ids
  )
  
  expect_equal(vary_result, 48)
  expect_equal(covy_result, 24)
})

test_that("pij preserves extremely small positive probabilities", {
  Ne_list <- list(c(3, 4), 3, c(1, 2), 1)
  p <- c(0.5, 0.5, 0.8, 1e-14)
  
  result <- pij(
    i = 1, j = 2,
    wi = 0, gi = 1,
    wj = 0, gj = 0,
    Ne_list = Ne_list, p = p
  )
  
  expected <- 0.5 * 0.5 * (1 - 0.8) * 1e-14
  
  expect_gt(result, 0)
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("Vary accounts for shared neighbors outside the leaf", {
  
  Ne_list <- list(3, 3, c(1, 2))
  Ne <- lengths(Ne_list)
  p <- rep(0.5, 3)
  
  W_full <- c(0, 0, 0)
  G_full <- c(0, 0, 0)
  Y_full <- c(2, 4, 6)
  
  unit_ids <- c(1, 2)
  
  result <- Vary(
    w = 0, g = 0,
    W = W_full[unit_ids],
    G = G_full[unit_ids],
    Y = Y_full[unit_ids],
    p = p,
    Ne = Ne,
    Ne_list = Ne_list,
    unit_ids = unit_ids
  )
  
  expect_equal(result, 368)
})


test_that("very small probabilities survive in marginal and joint exposure calculations", {
  p <- c(0.5, 0.5, 0.8, 1e-18)
  Ne_list <- list(c(3, 4), 3, c(1, 2), 1)

  # Marginal exposure with a single very-low-probability neighbor.
  p_marg <- c(0.5, 1e-18)
  Ne_marg <- list(2, 1)
  marg <- pi(i = 1, w = 0, g = 1, p = p_marg, Ne_list = Ne_marg)
  expect_gt(marg, 0)
  expect_equal(marg / ((1 - p_marg[1]) * p_marg[2]),
               1, tolerance = 1e-12)

  joint <- pij(
    i = 1, j = 2,
    wi = 0, gi = 1,
    wj = 0, gj = 0,
    Ne_list = Ne_list, p = p
  )
  expected <- (1 - p[1]) * (1 - p[2]) * (1 - p[3]) * p[4]
  expect_gt(joint, 0)
  expect_equal(joint / expected, 1, tolerance = 1e-12)

  # Nested neighborhoods make the reverse exposure combination impossible.
  zero_joint <- pij(
    i = 1, j = 2,
    wi = 0, gi = 0,
    wj = 0, gj = 1,
    Ne_list = Ne_list, p = p
  )
  expect_identical(zero_joint, 0)
})

test_that("NetworkCausalTree rejects malformed or directed adjacency matrices", {
  A <- matrix(0, nrow = 4, ncol = 4)
  A[1, 2] <- A[2, 1] <- 1
  A[3, 4] <- A[4, 3] <- 1

  run_tree <- function(adjacency) {
    NetworkCausalTree(
      X = matrix(c(0, 1, 0, 1), ncol = 1),
      Y = c(1, 2, 3, 4),
      W = c(0, 1, 0, 1),
      K = c(1, 1, 2, 2),
      A = adjacency,
      p = rep(0.5, 4),
      ratio_disc = 0.5,
      minsize = 1
    )
  }

  expect_error(run_tree(A[-1, , drop = FALSE]), "N x N")

  A_directed <- A
  A_directed[1, 2] <- 0
  expect_error(run_tree(A_directed), "symmetric")

  A_loop <- A
  A_loop[1, 1] <- 1
  expect_error(run_tree(A_loop), "zero diagonal")

  A_missing <- A
  A_missing[1, 2] <- NA_real_
  expect_error(run_tree(A_missing), "missing values")
})
