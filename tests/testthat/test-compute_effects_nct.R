test_that("compute_effects_nct runs on tiny tree and returns valid output", {
  set.seed(100)

  N <- 6
  W <- c(0,1,0,1,0,1)
  G <- c(0,0,1,1,0,1)
  Y <- 1 + W + G + rnorm(N, 0, 0.05)
  p <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)

  X <- matrix(1:6, ncol = 1)
  colnames(X) <- "X.1"

  nct_partition <- data.frame(
    NODE = c(1,2),
    OF = c(0,0.2),
    NOBS = c(N, 3),
    FILTER = c(NA, "data_tree$X.1 >= 3"),
    TERMINAL = c("PARENT", "LEAF")
  )

  output <- "estimation"
  minsize <- 1

  expect_warning(
    result <- compute_effects_nct(
      output = output,
      nct_partition = nct_partition,
      N = N,
      W = W, G = G, Y = Y, X = X,
      Ne = Ne, Ne_list = Ne_list,
      p = p, minsize = minsize
    ),
    "subpopulations not sufficiently represented"
  )

  expect_s3_class(result, "data.frame")

  expect_gte(nrow(result), 1)

  required_columns <- c("NODE","FILTER","TERMINAL","NOBS_TR","NOBS_EST",
                        "EFF1000_EST","SE1000_EST","CI1000_LOWER","CI1000_UPPER",
                        "EFF1101_EST","SE1101_EST","CI1101_LOWER","CI1101_UPPER",
                        "EFF1110_EST","SE1110_EST","CI1110_LOWER","CI1110_UPPER",
                        "EFF0100_EST","SE0100_EST","CI0100_LOWER","CI0100_UPPER")

  expect_true(all(required_columns %in% colnames(result)))

  expect_true(all(sapply(result$EFF1000_EST, is.numeric)))
  expect_true(all(sapply(result$SE1000_EST, is.numeric)))

  for (effect in c("1000", "1101", "1110", "0100")) {
    est <- result[[paste0("EFF", effect, "_EST")]]
    se <- result[[paste0("SE", effect, "_EST")]]

    lower <- result[[paste0("CI", effect, "_LOWER")]]
    upper <- result[[paste0("CI", effect, "_UPPER")]]

    expect_equal(lower, est - 1.96 * se)
    expect_equal(upper, est + 1.96 * se)
  }
})

