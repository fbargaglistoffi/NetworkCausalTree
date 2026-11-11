
test_that("K must be valid", {
  N <- 8
  K <- rep(1:2, each = 4)
  A <- matrix(0, N, N)
  X <- matrix(0, N, 1)
  Y <- rep(0, N)
  W <- rep(0, N)
  p <- rep(0.5, N)
  
  fit <- function(...) {
    args <- modifyList(
      list(X = X, Y = Y, W = W, A = A, K = K, p = p,
           ratio_disc = 0.5),
      list(...)
    )
    do.call(NetworkCausalTree, args)
  }
  
  expect_error(fit(K = NULL), "K must be")
  expect_error(fit(K = K[-1]), "K must be")
  expect_error(fit(K = replace(K, 1, NA)), "K must be")
})

test_that("A must not connect different clusters", {
  N <- 8
  K <- rep(1:2, each = 4)
  A <- matrix(0, N, N)
  A[1, 5] <- A[5, 1] <- 1
  
  expect_error(
    NetworkCausalTree(
      X = matrix(0, N, 1), Y = rep(0, N), W = rep(0, N),
      A = A, K = K, p = rep(0.5, N), ratio_disc = 0.5
    ),
    "block diagonal"
  )
})

test_that("p must contain valid probabilities", {
  N <- 8
  args <- list(
    X = matrix(0, N, 1), Y = rep(0, N), W = rep(0, N),
    A = matrix(0, N, N), K = rep(1:2, each = 4),
    ratio_disc = 0.5
  )
  
  for (bad_p in list(rep(0.5, N - 1), rep(0, N),
                     rep(1, N), rep(NA_real_, N),
                     rep(Inf, N))) {
    expect_error(
      do.call(NetworkCausalTree, c(args, list(p = bad_p))),
      "p must be"
    )
  }
})

test_that("method and output must be supported", {
  N <- 8
  args <- list(
    X = matrix(0, N, 1), Y = rep(0, N), W = rep(0, N),
    A = matrix(0, N, N), K = rep(1:2, each = 4),
    p = rep(0.5, N), ratio_disc = 0.5
  )
  
  expect_error(
    do.call(NetworkCausalTree, c(args, list(method = "invalid"))),
    "method must be"
  )
  expect_error(
    do.call(NetworkCausalTree, c(args, list(output = "invalid"))),
    "output must be"
  )
})

test_that("effect_weights must be nonnegative and valid", {
  N <- 8
  args <- list(
    X = matrix(0, N, 1), Y = rep(0, N), W = rep(0, N),
    A = matrix(0, N, N), K = rep(1:2, each = 4),
    p = rep(0.5, N), ratio_disc = 0.5
  )
  
  for (bad_weights in list(
    c(1.2, -0.2, 0, 0),
    c(1, 0, 0),
    c(NA_real_, 0, 0, 1),
    c(Inf, 0, 0, 0)
  )) {
    expect_error(
      do.call(NetworkCausalTree, c(args, list(effect_weights = bad_weights))),
      "effect_weights must"
    )
  }
})
