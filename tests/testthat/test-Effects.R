test_that("EffTau1000 computes correct effect on fake data", {
  N  <- 4
  W  <- c(1, 1, 0, 0)
  G  <- c(0, 1, 0, 1)
  Y  <- c(10, 20, 5, 7)
  p  <- rep(0.5, 4)
  Ne <- rep(1, 4)

  mock_pi <- function(...) rep(0.5, length(list(...)[[1]]))
  
  with_mocked_bindings(
    pi = mock_pi,
    {
      result <- EffTau1000(N, W, G, Y, p, Ne)
      expect_equal(result, 2.5, tolerance = 1e-8)
    }
  )
})

test_that("EffTau1101 computes correct effect on fake data", {
  N  <- 4
  W  <- c(1, 1, 0, 0)
  G  <- c(0, 1, 0, 1)
  Y  <- c(10, 20, 5, 7)
  p  <- rep(0.5, 4)
  Ne <- rep(1, 4)

  mock_pi <- function(...) rep(0.5, length(list(...)[[1]]))
  
  with_mocked_bindings(
    pi = mock_pi,
    {
      result <- EffTau1101(N, W, G, Y, p, Ne)
      expect_equal(result, 6.5, tolerance = 1e-8)
    }
  )
})

test_that("EffTau1110 computes correct effect on fake data", {
  N  <- 4
  W  <- c(1, 1, 0, 0)
  G  <- c(0, 1, 0, 1)
  Y  <- c(10, 20, 5, 7)
  p  <- rep(0.5, 4)
  Ne <- rep(1, 4)
  
  mock_pi <- function(...) rep(0.5, length(list(...)[[1]]))
  
  with_mocked_bindings(
    pi = mock_pi,
    {
      result <- EffTau1110(N, W, G, Y, p, Ne)
      expect_equal(result, 5, tolerance = 1e-8)
    }
  )
})

test_that("EffTau0100 computes correct effect on fake data", {
  N  <- 4
  W  <- c(1, 1, 0, 0)
  G  <- c(0, 1, 0, 1)
  Y  <- c(10, 20, 5, 7)
  p  <- rep(0.5, 4)
  Ne <- rep(1, 4)
  
  mock_pi <- function(...) rep(0.5, length(list(...)[[1]]))
  
  with_mocked_bindings(
    pi = mock_pi,
    {
      result <- EffTau0100(N, W, G, Y, p, Ne)
      expect_equal(result, 1, tolerance = 1e-8)
    }
  )
})

test_that("compute_population_effects returns correct vector on fake data", {
  N  <- 4
  W  <- c(1, 1, 0, 0)
  G  <- c(0, 1, 0, 1)
  Y  <- c(10, 20, 5, 7)
  p  <- rep(0.5, 4)
  Ne <- rep(1, 4)
  
  mock_pi <- function(...) rep(0.5, length(list(...)[[1]]))
  
  with_mocked_bindings(
    pi = mock_pi,
    {
      expected <- c(2.5, 6.5, 5, 1)
      
      result <- compute_population_effects(N, W, G, Y, p, Ne)
      
      expect_length(result, 4)
      expect_type(result, "double")
      expect_equal(result, expected, tolerance = 1e-8)
    }
  )
})

