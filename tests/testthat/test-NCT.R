test_that("plot_NCT runs without crashing on a small NCT object", {
  skip("Plot tests skipped (only tested interactively)")
  skip_if_not_installed("igraph")
  skip_if_not_installed("data.tree")

  NCT <- data.frame(
    FILTER = c("Root", "Root & X.1>=0.5"),
    NOBS_EST = c(10, 5),
    NOBS_TR = c(5, 2),
    EFF1000_EST = c(1.0, 2.0),
    SE1000_EST = c(0.1, 0.2),
    EFF1101_EST = c(0.5, 1.5),
    SE1101_EST = c(0.2, 0.3),
    EFF1110_EST = c(0.3, 1.0),
    SE1110_EST = c(0.1, 0.2),
    EFF0100_EST = c(0.2, 0.8),
    SE0100_EST = c(0.1, 0.1)
  )
  
  cov_names <- c("cov1")
  
  expect_silent(
    plot_NCT(NCT = NCT, cov_names = cov_names, title = "Test Tree")
  )
})
