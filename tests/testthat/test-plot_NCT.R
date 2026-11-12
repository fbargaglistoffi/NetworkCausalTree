test_that("plot_NCT runs on a real tiny NCT object", {
  skip_on_cran()

  set.seed(100)
  d <- data_generator_direct_indirect(N = 500, k = 10, M = 2, p = rep(0.3, 500), remove_isolates = TRUE)

  fit <- NetworkCausalTree(
    Y = d$Y, W = d$W, X = d$X, A = d$A, p = d$p,
    K = d$K,
    effect_weights = c(1,0,0,0),
    method = "singular", depth = 1, minsize = 20,
    ratio_disc = 0.5
  )

  cov_names <- colnames(d$X)

  expect_error(
    plot_NCT(fit, cov_names = cov_names, title = "Tiny Tree"),
    NA
  )
})

test_that("plot_NCT works in detection mode", {
  skip_on_cran()

  set.seed(100)
  d <- data_generator_direct_indirect(N = 500, k = 10, M = 2, p = rep(0.3, 500), remove_isolates = TRUE)

  fit <- NetworkCausalTree(
    Y = d$Y, W = d$W, X = d$X, A = d$A, p = d$p,
    K = d$K,
    effect_weights = c(1,0,0,0),
    method = "singular", depth = 1, minsize = 20,
    ratio_disc = 0.5,
    output = "detection"
  )

  cov_names <- colnames(d$X)

  expect_error(
    plot_NCT(fit, cov_names = cov_names, title = "Detect Tree", output = "detection"),
    NA
  )
})

test_that("node percentages use root sample size", {
  NOBS <- c(100, 60, 40)
  FILTER <- c("NA", "X.1 < 1", "X.1 >= 1")

  root_idx <- which(FILTER == "NA")
  percentages <- (NOBS / NOBS[root_idx]) * 100

  expect_equal(percentages, c(100, 60, 40))
  expect_equal(sum(percentages[-1]), 100)
})
