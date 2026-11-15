exposure_domain_fixture <- function() {
  N <- 16L
  neighbor <- as.vector(rbind(seq(2L, N, 2L), seq(1L, N, 2L)))
  A <- matrix(0, N, N)
  A[cbind(seq_len(N), neighbor)] <- 1
  W <- rep(c(0, 0, 1, 0, 0, 1, 1, 1), 2)
  G <- W[neighbor]
  list(X = matrix(0, N, 1), Y = 1 + 2 * W + 3 * G, W = W,
       A = A, K = rep(1:2, each = 8), p = rep(0.5, N),
       ratio_disc = 0.5, minsize = 1, depth = 1, output = "detection")
}

test_that("non-finite, negative, and weighted adjacency values are rejected", {
  for (bad in c(Inf, -Inf, -1, 0.5, 2)) {
    args <- exposure_domain_fixture()
    args$A[1, 2] <- args$A[2, 1] <- bad
    expect_error(do.call(NetworkCausalTree, args), "finite binary entries")
  }
  args <- exposure_domain_fixture()
  args$A[1, 2] <- args$A[2, 1] <- NA_real_
  expect_error(do.call(NetworkCausalTree, args), "missing values")
})

test_that("isolated units are rejected without changing the population", {
  for (output in c("detection", "estimation")) {
    for (all_isolated in c(FALSE, TRUE)) {
      args <- exposure_domain_fixture()
      args$output <- output
      if (all_isolated) {
        args$A[,] <- 0
      } else {
        args$A[1, 2] <- args$A[2, 1] <- 0
      }
      original <- args
      expect_error(do.call(NetworkCausalTree, args),
                   "isolated units.*zero probability.*no units have been removed")
      expect_identical(args, original)
    }
  }
})

test_that("supported dyads retain all four point estimates and uncertainty outputs", {
  for (output in c("detection", "estimation")) {
    args <- exposure_domain_fixture()
    args$output <- output
    set.seed(5)
    result <- do.call(NetworkCausalTree, args)
    expect_equal(result$NOBS_EST, 8)
    columns <- paste0("EFF", c("1000", "1101", "1110", "0100"), "_EST")
    estimates <- unlist(result[1, columns], use.names = FALSE)
    expect_equal(estimates, c(2, 2, 3, 3))
    if (output == "estimation") {
      expect_true(all(is.finite(as.matrix(result[, grep("^SE|^CI", names(result))]))))
    }
  }
})
