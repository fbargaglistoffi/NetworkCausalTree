test_that("Variance & covariance functions run on tiny network without errors", {
  set.seed(22)

  N <- 6
  W <- c(0,1,0,1,0,1)
  G <- c(0,0,1,1,0,1)
  Y <- 1 + W + G + rnorm(N, 0, 0.05)
  p <- rep(0.5, N)
  Ne <- rep(1, N)
  Ne_list <- list(2,1,4,3,6,5)

  pi  <- NetworkCausalTree:::pi
  pij <- NetworkCausalTree:::pij
  shared_neigh <- NetworkCausalTree:::shared_neigh

  expect_true({
    v00 <- Vary(N, w=0, g=0, Y, W, G, p, Ne, Ne_list)
    is.numeric(v00) || is.na(v00)
  })

  expect_true({
    c0010 <- Covy(0,0,1,0, N, Y, W, G, p, Ne, Ne_list)
    is.numeric(c0010) || is.na(c0010)
  })

  expect_true(is.numeric(Vartau1000(N,Y,W,G,p,Ne,Ne_list)) || is.na(Vartau1000(N,Y,W,G,p,Ne,Ne_list)))
  expect_true(is.numeric(Vartau1101(N,Y,W,G,p,Ne,Ne_list)) || is.na(Vartau1101(N,Y,W,G,p,Ne,Ne_list)))
  expect_true(is.numeric(Vartau1110(N,Y,W,G,p,Ne,Ne_list)) || is.na(Vartau1110(N,Y,W,G,p,Ne,Ne_list)))
  expect_true(is.numeric(Vartau0100(N,Y,W,G,p,Ne,Ne_list)) || is.na(Vartau0100(N,Y,W,G,p,Ne,Ne_list)))
})

