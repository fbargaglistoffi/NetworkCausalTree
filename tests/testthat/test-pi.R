test_that("pi computes correct marginal exposure probabilities", {
  p  <- c(0.5, 0.8, 0.2)
  Ne <- c(2, 1, 3)
  
  # spot check
  # confirm pi() works correctly
  expected_1 <- 0.5 * (1 - 0.5)^2
  result_1 <- pi(i = 1, w = 1, g = 0, p = p, Ne = Ne)
  expect_equal(result_1, expected_1, tolerance = 1e-8)
  
  expected_2 <- (1 - 0.8) * (1 - (1 - 0.8)^1)  
  result_2 <- pi(i = 2, w = 0, g = 1, p = p, Ne = Ne)
  expect_equal(result_2, expected_2, tolerance = 1e-8)
  
  # batch check
  # makes sure it works correclty and consistent for multiple i values
  results <- sapply(1:3, function(i) pi(i, w = 1, g = 0, p = p, Ne = Ne))
  expected <- c(0.5*(0.5^2), 0.8*(0.2^1), 0.2*(0.8^3))
  expect_equal(results, expected, tolerance = 1e-8)
})
