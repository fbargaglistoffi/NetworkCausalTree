test_that("generate_clustered_networks produces valid adjacency matrices", {
  set.seed(123)
  
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
})

test_that("expand.grid.unique works correctly", {
  x <- c(1, 2, 3)
  y <- c(2, 3, 4)
  
  result <- expand.grid.unique(x, y, include.equals = FALSE)
  
  # check structure
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 2)
  
  # ensure no equal elements in rows
  expect_true(all(result[,1] != result[,2]))
  
  # include.equals = TRUE allows equal values
  result2 <- expand.grid.unique(x, y, include.equals = TRUE)
  expect_true(any(result2[,1] == result2[,2]))
})

test_that("shared_neigh computes shared neighbors accurately", {
  # Example adjacency and neighbor list
  A <- matrix(c(
    0,1,1,
    1,0,1,
    1,1,0
  ), nrow = 3, byrow = TRUE)
  
  Ne_list <- list(
    c(2,3),
    c(1,3),
    c(1,2)
  )
  
  expect_equal(shared_neigh(1, 2, Ne_list), 1)  # shared neighbor = node 3
  expect_equal(shared_neigh(1, 3, Ne_list), 1)  # shared = node 2
  expect_equal(shared_neigh(2, 3, Ne_list), 1)  # shared = node 1
})
