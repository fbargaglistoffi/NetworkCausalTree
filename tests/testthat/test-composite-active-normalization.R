# Four disjoint dyads provide two units in every observed exposure cell.
# With p = 1/2 and Y = 2 W + 3 G, the four HT contrasts are (2, 2, 3, 3).
composite_normalization_fixture <- function() {
  neighbors <- as.list(c(2L, 1L, 4L, 3L, 6L, 5L, 8L, 7L))
  W <- c(0, 0, 1, 0, 0, 1, 1, 1)
  G <- W[unlist(neighbors)]
  list(method = "composite", alpha = 0.5, beta = 0.5,
       gamma = 0, delta = 0, N = 8L, W = W, G = G,
       Y = 2 * W + 3 * G, p = rep(0.5, 8), Ne = rep(1, 8),
       Ne_list = neighbors, total_variance = NULL, nleafs = 1)
}

test_that("inactive zero and non-finite normalizers do not affect the criterion", {
  args <- composite_normalization_fixture()
  for (unused in list(c(0, 0), c(NA_real_, Inf), c(NaN, -Inf))) {
    args$population_effects <- c(1, 2, unused)
    expect_equal(do.call(NetworkCausalTree:::compute_OF_Value, args), 2.5)
  }
})

test_that("selected zero and non-finite normalizers still produce an error", {
  args <- composite_normalization_fixture()
  args$population_effects <- c(0, 2, 1, 1)
  expect_error(do.call(NetworkCausalTree:::compute_OF_Value, args),
               "zero for selected effects")
  for (bad in c(NA_real_, Inf, NaN)) {
    args$population_effects <- c(1, bad, 1, 1)
    expect_error(do.call(NetworkCausalTree:::compute_OF_Value, args),
                 "non-finite values for selected effects")
  }
})

test_that("all-active composite scores retain the weighted normalized formula", {
  args <- composite_normalization_fixture()
  args$alpha <- args$beta <- args$gamma <- args$delta <- 0.25
  args$population_effects <- c(1, 2, -3, 6)
  expected <- sum(0.25 * (c(2, 2, 3, 3) / args$population_effects)^2)
  expect_equal(do.call(NetworkCausalTree:::compute_OF_Value, args), expected)
})
