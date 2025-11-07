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
  
  # ERGM network is available but not covered in the test suite because of
  # computational reasons (as it uses the ergm package internally)
})

test_that("expand.grid.unique works correctly", {
  x <- c(1, 2, 3)
  y <- c(2, 3, 4)
  
  result <- expand.grid.unique(x, y, include.equals = FALSE)
  
  # result must be a dataframe and have only 2 columns 
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 2)
  
  # because include.equals = FALSE, there will be no rows with equal elements
  expect_true(all(result[,1] != result[,2]))
  
  # include.equals = TRUE allows equal values
  result2 <- expand.grid.unique(x, y, include.equals = TRUE)
  expect_true(any(result2[,1] == result2[,2]))
})

test_that("shared_neigh counts common neighbors correctly", {
  
  Ne_list <- list(
    c(2),
    c(1,3),
    c(2)
  )
  
  expect_equal(shared_neigh(1, 2, Ne_list), 0)
  
  expect_equal(shared_neigh(2, 3, Ne_list), 0)
  
  expect_equal(shared_neigh(1, 3, Ne_list), 1)
  
  expect_equal(shared_neigh(2, 2, Ne_list), 2)
})
