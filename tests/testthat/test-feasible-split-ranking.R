feasible_split_fixture <- function(method = "singular") {
  N <- 32L
  W <- rep(c(0, 0, 1, 0, 0, 1, 1, 1), 4)
  neighbors <- as.list(as.vector(rbind(seq(2L, N, 2L), seq(1L, N, 2L))))
  G <- W[unlist(neighbors)]
  X <- matrix(rep(1:4, each = 8), ncol = 1,
              dimnames = list(NULL, "X.1"))
  tau <- c(10, -2, -2, -2)[X[, 1]]
  list(method = method,
       alpha = if (method == "composite") 0.5 else 1,
       beta = if (method == "composite") 0.5 else 0,
       gamma = 0, delta = 0, W = W, G = G, Y = tau * W, X = X,
       p = rep(0.5, N), Ne = rep(1L, N), Ne_list = neighbors,
       population_effects = rep(1, 4), total_variance = 0,
       nleafs = 1, N_tr = N, N_est = N)
}

test_that("a lower-scoring feasible split is considered before stopping", {
  for (method in c("singular", "composite")) {
    args <- feasible_split_fixture(method)
    loose <- do.call(NetworkCausalTree:::compute_OF_Split, args)
    expect_equal(as.numeric(loose["split"]), 2)
    expect_equal(as.numeric(loose["of"]), 28)
    args$minsize <- 3
    feasible <- do.call(NetworkCausalTree:::compute_OF_Split, args)
    expect_equal(as.numeric(feasible["split"]), 3)
    expect_equal(as.numeric(feasible["of"]), 10)

    # Exercise the actual caller, not only the scoring helper.
    args$total_variance <- args$nleafs <- NULL
    args$N <- 32L
    args$depth <- 1
    tree <- do.call(NetworkCausalTree:::identify_partitions_nct, args)
    expect_equal(nrow(tree), 3L)
    expect_equal(tree$NOBS[-1], c(16L, 16L))
    expect_setequal(tree$FILTER[-1], c("X.1 >= 3", "X.1 < 3"))
  }
})

test_that("all criteria screen both children and accept the minsize boundary", {
  for (method in c("singular", "composite", "penalized")) {
    args <- feasible_split_fixture(method)
    args$minsize <- 4
    result <- do.call(NetworkCausalTree:::compute_OF_Split, args)
    expect_equal(as.numeric(result["split"]), 3)
    args$minsize <- 5
    result <- do.call(NetworkCausalTree:::compute_OF_Split, args)
    expect_true(all(is.na(result)))
  }
})
