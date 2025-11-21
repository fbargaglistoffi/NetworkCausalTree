test_that("generate_clustered_networks produces valid adjacency matrices", {
  set.seed(100)
  
  N <- 20
  k <- 4
  M <- 1
  X <- matrix(rbinom(N, 1, 0.5), ncol = M)
  
  # ER network
  A_er <- generate_clustered_networks(
    k = k, N = N,
    method_networks = "er",
    param_er = 0.2,
    var_homophily_ergm = NULL,
    coef_ergm = NULL,
    X = X
  )
  
  expect_true(is.matrix(A_er))
  expect_equal(dim(A_er), c(N, N))
  expect_true(all(A_er %in% c(0, 1)))
  expect_equal(sum(diag(A_er)), 0) # no self loops
  
  # Barabasi network
  A_sf <- generate_clustered_networks(
    k = k, N = N,
    method_networks = "sf",
    param_er = 0.2,
    var_homophily_ergm = NULL,
    coef_ergm = NULL,
    X = X
  )
  
  expect_true(is.matrix(A_sf))
  expect_equal(dim(A_sf), c(N, N))
  expect_true(all(A_sf %in% c(0, 1)))
  
  # ERGM network is available but not covered in the test suite because of
  # computational reasons (as it uses the ergm package internally)
})

