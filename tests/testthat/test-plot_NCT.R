test_that("plot_NCT runs on a real tiny NCT object", {
  skip_on_cran()
  
  set.seed(100)
  d <- data_generator(N = 100, k = 5, M = 2, p = rep(0.3,100), remove_isolates = TRUE)
  
  fit <- NetworkCausalTree(
    Y = d$Y, W = d$W, X = d$X, A = d$A, p = d$p,
    method = "OF", depth = 1, minsize = 20,
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
  d <- data_generator(N = 100, k = 5, M = 2, p = rep(0.3,100), remove_isolates = TRUE)
  
  fit <- NetworkCausalTree(
    Y = d$Y, W = d$W, X = d$X, A = d$A, p = d$p,
    method = "OF", depth = 1, minsize = 20,
    ratio_disc = 0.5,
    output = "detection"
  )
  
  cov_names <- colnames(d$X)
  
  expect_error(
    plot_NCT(fit, cov_names = cov_names, title = "Detect Tree", output = "detection"),
    NA
  )
})
