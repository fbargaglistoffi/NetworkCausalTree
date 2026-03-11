test_that("pij computes correct joint probabilities for no treatment, no spillover", {
  
  Ne_list <- list(
    c(2),
    c(1,3),
    c(2)
  )
  
  Ne <- c(1,2,1)
  p  <- c(0.5, 0.8, 0.2)
  
  expected12 <- (1 - p[1])^(1 + Ne[1]) * (1 - p[2])^(1 + Ne[2]) # ... = 0.002
  result12   <- pij(1, 2, 0, 0, 0, 0, Ne, Ne_list, p)
  expect_equal(result12, expected12, tolerance = 1e-8)
  
  expected23 <- (1 - p[2])^(1 + Ne[2]) * (1 - p[3])^(1 + Ne[3]) # ... = 0.0256
  result23   <- pij(2, 3, 0, 0, 0, 0, Ne, Ne_list, p)
  expect_equal(result23, expected23, tolerance = 1e-8)
})

test_that("pij respects parameter order with distinct wi/gi/wj/gj", {
  Ne_list <- list(c(2), c(1))
  Ne <- c(2, 1)
  p <- c(0.3, 0.7)
  
  result1 <- pij(1, 2, wi = 1, gi = 0, wj = 0, gj = 1, Ne = Ne, Ne_list = Ne_list, p = p)
  result2 <- pij(1, 2, wi = 0, gi = 1, wj = 1, gj = 0, Ne = Ne, Ne_list = Ne_list, p = p)
  
  expect_true(is.numeric(result1))
  expect_true(is.numeric(result2))
  expect_false(isTRUE(all.equal(result1, result2)))
})
