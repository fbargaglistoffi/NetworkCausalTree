test_that("Covy runs on tiny network without errors", {
  set.seed(100)
  
  N  <- 6
  W  <- c(0,1,0,1,0,1)
  G  <- c(0,0,1,1,0,1)
  Y  <- 1 + W + G + rnorm(N, 0, 0.05)
  p  <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)
  
  c0010 <- Covy(0,0,1,0, N, Y, W, G, p, Ne, Ne_list)
  
  expect_true(is.numeric(c0010) || is.na(c0010))
})