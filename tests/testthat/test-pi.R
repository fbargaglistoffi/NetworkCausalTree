test_that("pi computes correct marginal exposure probabilities", {
  p <- c(0.5, 0.8, 0.2)
  
  Ne_list <- list(
    c(2, 3),
    c(1),
    c(1, 2)
  )
  
  expected_1 <- 0.5 * (1 - 0.8) * (1 - 0.2)
  result_1 <- pi(i = 1, w = 1, g = 0, p = p, Ne_list = Ne_list)
  expect_equal(result_1, expected_1, tolerance = 1e-8)
  
  expected_2 <- (1 - 0.8) * (1 - (1 - 0.5))
  result_2 <- pi(i = 2, w = 0, g = 1, p = p, Ne_list = Ne_list)
  expect_equal(result_2, expected_2, tolerance = 1e-8)
  
  results <- pi(i = 1:3, w = 1, g = 0, p = p, Ne_list = Ne_list)
  
  expected <- c(
    0.5 * (1 - 0.8) * (1 - 0.2),
    0.8 * (1 - 0.5),
    0.2 * (1 - 0.5) * (1 - 0.8)
  )
  
  expect_equal(results, expected, tolerance = 1e-8)
})

test_that("pi marginal exposure probabilities sum to one", {
  Ne_list <- list(
    c(2, 3),
    c(1, 3),
    c(1, 2)
  )
  
  p <- c(0.3, 0.7, 0.4)
  
  for (i in 1:3) {
    total <- sum(
      sapply(0:1, function(w) {
        sum(sapply(0:1, function(g) {
          pi(i, w, g, p, Ne_list)
        }))
      })
    )
    
    expect_equal(total, 1, tolerance = 1e-10)
  }
})
