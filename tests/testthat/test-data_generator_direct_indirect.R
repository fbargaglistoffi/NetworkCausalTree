test_that("data_generator_direct_indirect runs and returns valid structure", {
  set.seed(100)
  
  N <- 50
  M <- 3
  k <- 5
  p <- rep(0.3, N)
  
  data <- data_generator_direct_indirect(
    N = N, M = M, k = k, p = p,
    method_networks = "er", param_er = 0.2
  )
  
  expect_type(data, "list")
  
  expect_named(
    data,
    c("X","Y","W","A","G","K","p"),
    ignore.order = TRUE
  )
  
  n_new <- length(data$W)
  expect_gt(n_new, 0)
  expect_equal(n_new, nrow(data$X))
  expect_equal(n_new, length(data$Y))
  expect_equal(n_new, length(data$G))
  expect_equal(n_new, length(data$K))
  expect_equal(n_new, length(data$p))
  expect_equal(n_new, nrow(data$A))
  expect_equal(n_new, ncol(data$A))
  
  expect_equal(ncol(data$X), M)
  
  expect_true(all(data$W %in% c(0,1)))
  expect_true(all(data$G %in% c(0,1)))
  
  expect_true(is.matrix(data$A))
  expect_true(all(data$A %in% c(0,1)))
  expect_equal(sum(diag(data$A)), 0)
  
  # remove_isolates = TRUE
  expect_true(all(rowSums(data$A) > 0))
})

test_that("data_generator_direct_indirect handles homogeneous vs heterogeneous effects", {
  set.seed(100)
  N <- 80
  
  p <- rep(0.5, N)
  
  data_het <- data_generator_direct_indirect(N = N, p = p, het = TRUE)
  data_hom <- data_generator_direct_indirect(N = N, p = p, het = FALSE)
  
  expect_type(data_het, "list")
  expect_type(data_hom, "list")
  
  sd_het <- sd(data_het$Y)
  sd_hom <- sd(data_hom$Y)
  
  expect_gt(abs(sd_het - sd_hom), 0.1)
})

test_that("data_generator_direct_indirect errors if p length mismatches N", {
  expect_error(data_generator_direct_indirect(N = 20, p = rep(0.2, 10)))
})

test_that("data_generator_direct_indirect works with non-default N and no p", {
  result <- data_generator_direct_indirect(N = 100, k = 10, remove_isolates = FALSE)
  expect_equal(length(result$p), 100)
})

test_that("data_generator_direct_indirect errors when N not divisible by k", {
  expect_error(data_generator_direct_indirect(N = 101, k = 10))
})
