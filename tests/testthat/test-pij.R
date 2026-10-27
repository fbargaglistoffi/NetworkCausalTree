test_that("pij computes correct joint probabilities for no treatment, no spillover", {
  
  Ne_list <- list(
    c(2),
    c(1, 3),
    c(2)
  )
  
  p <- c(0.5, 0.8, 0.2)
  
  expected12 <- (1 - p[1]) * (1 - p[2]) * (1 - p[3])
  result12 <- pij(1, 2, 0, 0, 0, 0, Ne_list, p)
  expect_equal(result12, expected12, tolerance = 1e-8)
  
  expected23 <- (1 - p[1]) * (1 - p[2]) * (1 - p[3])
  result23 <- pij(2, 3, 0, 0, 0, 0, Ne_list, p)
  expect_equal(result23, expected23, tolerance = 1e-8)
})

test_that("pij respects parameter order with distinct wi/gi/wj/gj", {
  Ne_list <- list(c(2), c(1))
  p <- c(0.3, 0.7)
  
  result1 <- pij(1, 2, wi = 1, gi = 0, wj = 0, gj = 1,
                 Ne_list = Ne_list, p = p)
  
  result2 <- pij(1, 2, wi = 0, gi = 1, wj = 1, gj = 0,
                 Ne_list = Ne_list, p = p)
  
  expect_equal(result1, 0.09, tolerance = 1e-8)
  expect_equal(result2, 0.49, tolerance = 1e-8)
})

test_that("pij handles shared neighbors correctly", {
  Ne_list <- list(
    c(3),
    c(3),
    c(1, 2)
  )
  
  p <- rep(0.5, 3)
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 1,
    wj = 1, gj = 1,
    Ne_list = Ne_list, p = p
  )
  
  expect_equal(result, 0.125, tolerance = 1e-8)
})

test_that("pij returns zero for impossible exposure combinations", {
  Ne_list <- list(
    c(2),
    c(1)
  )
  
  p <- rep(0.5, 2)
  
  result <- pij(
    i = 1, j = 2,
    wi = 1, gi = 0,
    wj = 0, gj = 0,
    Ne_list = Ne_list, p = p
  )
  
  expect_equal(result, 0, tolerance = 1e-8)
})

test_that("pij matches exact enumeration with heterogeneous probabilities", {
  
  Ne_list <- list(
    c(3),
    c(3),
    c(1, 2)
  )
  
  p <- c(0.3, 0.7, 0.4)
  
  assignments <- as.matrix(expand.grid(rep(list(0:1), 3)))
  
  probs <- apply(assignments, 1, function(W) {
    prod(p^W * (1 - p)^(1 - W))
  })
  
  G1 <- apply(assignments, 1, function(W) {
    as.integer(any(W[Ne_list[[1]]] == 1))
  })
  
  G2 <- apply(assignments, 1, function(W) {
    as.integer(any(W[Ne_list[[2]]] == 1))
  })
  
  matches <- assignments[, 1] == 1 &
    assignments[, 2] == 1 &
    G1 == 1 & G2 == 1
  
  exact <- sum(probs[matches])
  
  analytical <- pij(
    i = 1, j = 2,
    wi = 1, gi = 1,
    wj = 1, gj = 1,
    Ne_list = Ne_list, p = p
  )
  
  expect_equal(analytical, exact, tolerance = 1e-10)
  expect_equal(exact, 0.084, tolerance = 1e-10)
})

test_that("pij matches enumeration for all 16 exposure combinations", {
  Ne_list <- list(c(2, 3), c(1, 3), c(1, 2))
  p <- c(0.3, 0.7, 0.4)
  
  assignments <- as.matrix(expand.grid(rep(list(0:1), 3)))
  
  probs <- apply(assignments, 1, function(W) {
    prod(p^W * (1 - p)^(1 - W))
  })
  
  G1 <- apply(assignments, 1, function(W) {
    as.integer(any(W[Ne_list[[1]]] == 1))
  })
  
  G2 <- apply(assignments, 1, function(W) {
    as.integer(any(W[Ne_list[[2]]] == 1))
  })
  
  for (wi in 0:1) {
    for (gi in 0:1) {
      for (wj in 0:1) {
        for (gj in 0:1) {
          
          matches <- assignments[, 1] == wi &
            assignments[, 2] == wj &
            G1 == gi & G2 == gj
          
          exact <- sum(probs[matches])
          
          analytical <- pij(
            1, 2, wi, gi, wj, gj,
            Ne_list = Ne_list, p = p
          )
          
          expect_equal(analytical, exact, tolerance = 1e-10)
        }
      }
    }
  }
})
