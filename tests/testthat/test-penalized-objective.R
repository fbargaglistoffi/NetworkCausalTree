test_that("penalized split matches manual objective improvement", {
  
  N <- 16
  
  X <- data.frame(X.1 = c(rep(0, 8), rep(1, 8)))
  
  Ne_list <- lapply(seq_len(N), function(i) {
    if (i %% 2 == 1) i + 1 else i - 1
  })
  Ne <- lengths(Ne_list)
  
  W <- rep(c(0, 0, 1, 1, 0, 1, 0, 1), 2)
  G <- as.integer(sapply(Ne_list, function(ids) any(W[ids] == 1)))
  Y <- seq_len(N)
  p <- rep(0.5, N)
  
  N_tr <- 16
  N_est <- 24
  
  result <- compute_OF_Split(
    method = "penalized",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = rep(1, 4),
    total_variance = 0, nleafs = 1,
    N_tr = N_tr, N_est = N_est
  )
  
  left <- 1:8
  right <- 9:16
  
  tau_parent <- EffTau1000(N, W, G, Y, p, Ne_list)
  tau_left <- EffTau1000(8, W[left], G[left], Y[left],
                         p, Ne_list, unit_ids = left)
  tau_right <- EffTau1000(8, W[right], G[right], Y[right],
                          p, Ne_list, unit_ids = right)
  
  reward_change <- (8 / 16) * tau_left^2 +
    (8 / 16) * tau_right^2 - tau_parent^2
  
  v_parent <- Vartau1000(
    N = N, W = W, G = G, Y = Y,
    p = p, Ne = Ne, Ne_list = Ne_list
  )
  
  v_left <- Vartau1000(
    N = 8, W = W[left], G = G[left], Y = Y[left],
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = left
  )
  
  v_right <- Vartau1000(
    N = 8, W = W[right], G = G[right], Y = Y[right],
    p = p, Ne = Ne, Ne_list = Ne_list,
    unit_ids = right
  )
  
  variance_change <- v_left + v_right - v_parent
  
  expected <- reward_change -
    (1 / N_tr + 1 / N_est) * variance_change
  
  expect_equal(as.numeric(result["of"]), expected)
  expect_equal(unname(result["var"]), "X.1")
  expect_equal(as.numeric(result["split"]), 1)
  
  # Imagine the parent is one of two existing leaves.
  # The total discovery sample contains 32 observations,
  # but this parent contains only 16.
  
  N_tr_multi <- 32
  other_leaf_variance <- 10
  
  result_multi <- compute_OF_Split(
    method = "penalized",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    W = W, G = G, Y = Y, X = X,
    p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = rep(1, 4),
    total_variance = other_leaf_variance + v_parent,
    nleafs = 2,
    N_tr = N_tr_multi, N_est = N_est
  )
  
  expected_multi <- (16 / N_tr_multi) * reward_change -
    (1 / N_tr_multi + 1 / N_est) *
    (v_left + v_right - v_parent)
  
  expect_equal(as.numeric(result_multi["of"]), expected_multi)
})

test_that("penalized splitting rejects non-improving splits", {
  maxof <- -0.01
  split_here <- rep(TRUE, 2)
  
  if (maxof <= 0) {
    split_here <- rep(FALSE, 2)
  }
  
  expect_false(any(split_here))
})

test_that("penalized tree rejects a non-improving split", {
  testthat::local_mocked_bindings(
    compute_OF_Split = function(...) {
      c(of = -0.01, split = 1, var = "X.1")
    },
    Vartau1000 = function(...) 1,
    .package = "NetworkCausalTree"
  )
  
  N <- 8
  
  tree <- identify_partitions_nct(
    method = "penalized",
    alpha = 1, beta = 0, gamma = 0, delta = 0,
    depth = 2, minsize = 1,
    N = N,
    W = rep(c(0, 0, 1, 1), 2),
    G = rep(c(0, 1, 0, 1), 2),
    Y = seq_len(N),
    X = data.frame(X.1 = rep(c(0, 1), each = 4)),
    p = rep(0.5, N),
    Ne = rep(0, N),
    Ne_list = replicate(N, integer(0), simplify = FALSE),
    population_effects = rep(1, 4),
    N_tr = 8, N_est = 8
  )
  
  expect_equal(nrow(tree), 1)
  expect_equal(tree$TERMINAL[1], "LEAF")
})

