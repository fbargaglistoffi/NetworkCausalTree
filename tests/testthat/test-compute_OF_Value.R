test_that("compute_OF_Value works for singular case on tiny example", {

  Ne_list <- list(
    c(2),
    c(1,3),
    c(2)
  )
  
  Ne <- c(1,2,1)
  p  <- c(0.5, 0.8, 0.2)

  W <- c(1,0,1)
  G <- c(0,1,0)

  Y <- c(10, 5, 15)
  
  N <- length(Y)

  # not actually used but need for argument
  pe <- compute_population_effects(N, W, G, Y, p, Ne)

  alpha <- 1
  beta <- 0
  gamma <- 0
  delta <- 0
  
  result <- compute_OF_Value(
    method = "singular",
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    N = N, W = W, G = G, Y = Y, p = p, Ne = Ne, Ne_list = Ne_list,
    population_effects = pe,
    total_variance = NULL, nleafs = NULL
  )
  
  expected <- (EffTau1000(N, W, G, Y, p, Ne))^2
  
  expect_equal(result, expected, tolerance = 1e-8)
})

test_that("compute_OF_Value stops when population_effects contains zero in composite method", {
  expect_error(
    compute_OF_Value("composite", 1, 0, 0, 0, N = 3,
                     W = c(0,1,0), G = c(0,0,1), Y = c(1,2,3),
                     p = rep(0.5,3), Ne = rep(1,3), Ne_list = NULL,
                     population_effects = c(0, 1, 1, 1),
                     total_variance = NULL, nleafs = 1), "population_effects contains zero"
  )
})

test_that("compute_OF_Value stops on unknown method", {
  expect_error(
    compute_OF_Value("unknown", 1, 0, 0, 0, N = 3,
                     W = c(0,1,0), G = c(0,0,1), Y = c(1,2,3),
                     p = rep(0.5,3), Ne = rep(1,3), Ne_list = NULL,
                     population_effects = c(1,1,1,1),
                     total_variance = NULL, nleafs = 1), "Unknown method"
  )
})

test_that("compute_OF_Value stops when population_effects contains NA in composite method", {
  expect_error(
    compute_OF_Value("composite", 1, 0, 0, 0, N = 3,
                     W = c(0,1,0), G = c(0,0,1), Y = c(1,2,3),
                     p = rep(0.5,3), Ne = rep(1,3), Ne_list = NULL,
                     population_effects = c(NA, 1, 1, 1),
                     total_variance = NULL, nleafs = 1), "non-finite values"
  )
})

test_that("compute_OF_Value stops when population_effects is wrong length in composite method", {
  expect_error(
    compute_OF_Value("composite", 1, 0, 0, 0, N = 3,
                     W = c(0,1,0), G = c(0,0,1), Y = c(1,2,3),
                     p = rep(0.5,3), Ne = rep(1,3), Ne_list = NULL,
                     population_effects = c(1, 1, 1),
                     total_variance = NULL, nleafs = 1), "length 4"
  )
})
