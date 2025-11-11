test_that("NetworkCausalTree runs and returns expected structure", {
  set.seed(100)
  N <- 20
  M <- 2
  
  dataset_list <- list(
    X = matrix(rbinom(N, size = 1, prob = 0.5)),
    Y = rnorm(N),
    W = rbinom(N, 1, 0.5),
    A = {
      A <- matrix(sample(0:1, N * N, replace = TRUE, prob = c(0.8, 0.2)), nrow = N)
      diag(A) <- 0
      A
    },
    K = sample(1:4, N, replace = TRUE),
    p = rbinom(N, size = 1, prob = 0.5)
  )
  
  result <- NetworkCausalTree(
    X = dataset_list[["X"]],
    Y = dataset_list[["Y"]],
    W = dataset_list[["W"]],
    A = dataset_list[["A"]],
    K = dataset_list[["K"]],
    p = dataset_list[["p"]],
    effect_weights = c(1, 0, 0, 0),
    ratio_disc = 0.5,
    depth = 3,
    minsize = 5,
    method = "singular",
    output = "estimation"
  )
  
  expect_true(is.data.frame(result))
  expect_true(all(c("FILTER", "NOBS_EST") %in% names(result)))
  expect_gt(nrow(result), 0)
})
