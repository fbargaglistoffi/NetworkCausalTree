test_that("caller covariate names do not change discovery or held-out estimation", {
  N <- 64L
  K <- rep(1:4, each = 16)
  neighbors <- as.vector(rbind(seq(2L, N, 2L), seq(1L, N, 2L)))
  A <- matrix(0, N, N)
  A[cbind(seq_len(N), neighbors)] <- 1
  W <- rep(c(0, 0, 1, 0, 0, 1, 1, 1), 8)
  G <- W[neighbors]
  X <- cbind(rep(rep(c(0, 1), each = 8), 4), rep(0, N))
  colnames(X) <- NULL
  Y <- 2 + (2 + 6 * X[, 1]) * W + G
  fit <- function(x, output) {
    set.seed(102)
    NetworkCausalTree(X = x, Y = Y, W = W, A = A, K = K,
                     p = rep(0.5, N), ratio_disc = 0.5,
                     depth = 1, minsize = 2, output = output)
  }
  labels <- list(c("X.age", "X.bmi"), c("age", "X.bmi"),
                 c("X.age", "X.age"), c("X.blood pressure", "X.risk-score"),
                 c("X.2", "X.1"))
  for (output in c("detection", "estimation")) {
    reference <- fit(X, output)
    expect_equal(nrow(reference), 3L)
    expect_equal(reference$NOBS_EST, c(32, 16, 16))
    expect_setequal(reference$FILTER[-1], c("X.1 >= 1", "X.1 < 1"))
    for (names_in in labels) {
      named <- X
      colnames(named) <- names_in
      expect_equal(fit(named, output), reference)
      # Internal renaming must not mutate labels retained by the caller for plotting.
      expect_identical(colnames(named), names_in)
      frame <- as.data.frame(named)
      expect_equal(fit(frame, output), reference)
    }
  }
})
