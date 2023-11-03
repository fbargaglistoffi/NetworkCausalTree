#' @title
#' Synthetic data generator
#'
#' @description
#' Generates Network Causal Tree synthetic data.
#'
#' @param  N Sample size (default: 2000).
#' @param  K Number of binary regressors (default: 5).
#' @param  m Number of clusters (default: 40).
#' @param  p  N x 1 vector, Probability to be assigned to the active individual
#' intervention (default: rep(0.2,2000))
#' @param het TRUE if the treatment effects 1000 and 1101 are heterogeneous with
#' respect to the first regressor (+taui with X1=0, -taui with X0=0), FALSE if
#' constant (+taui) (default: TRUE).
#' @param h Absolute value of the treatment effects 1000 and 1101
#' (default: 2).
#' @param  method_networks Method to generate the m networks:
#' "ergm" (Exponential Random Graph Models), "er" (Erdos Renyi), "sf"
#' (Barabasi-Albert model) (default: "er").
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er Probability of the "er" model, if used (default: 0.2).
#' @param  coef_ergm Coefficients of the "ergm" model , if used (default: NULL).
#' @param  var_homophily_ergm Variable to account for homophily in the "ergm"
#' model (default: NULL).
#'
#' @return A list of synthetic data containing:
#' - NxK covariates matrix (`X`).
#' - Nx1 outcome vector (`Y`),
#' - Nx1 individual intervention vector (`W`),
#' - NxN adjacency matrix (`A`),
#' - Nx1 neighborhood intervention vector (`G`),
#' - Nx1 group membership vector (`M`),
#' - Nx1 probability to be assigned to the active individual intervention vector
#' (`p`),
#'
#' @import igraph
#'
#' @export

data_generator = function(N = 2000,
                          K = 5,
                          m = 40,
                          p = rep(0.2,2000),
                          het = TRUE,
                          h = 2,
                          method_networks = "er",
                          param_er = 0.1,
                          coef_ergm = NULL,
                          var_homophily_ergm = NULL){

  # Generate Covariates
  X <- NULL
  for (k in 1 : K) {
    x <- rbinom(N, 1, 0.5)
    X <- cbind(X, x)
    colnames(X)[k] <- paste0(colnames(X)[k], k)
  }

  # Generate m networks
  A <- generate_clustered_networks(N = N,
                                   m = m,
                                   method_networks = method_networks,
                                   param_er = param_er,
                                   coef_ergm = coef_ergm,
                                   var_homophily_ergm = var_homophily_ergm,
                                   X = X)
  
  net <- igraph::graph_from_adjacency_matrix(A)

  # Group Indicator
  cluster_size <- N / m
  M <- c(rep(1 : m, cluster_size))
  M <- sort(M)
  levels(M) <- c(1 : m)
  

  # Randomly assign unit to treatment arms
  W <- rbinom(N, 1, prob = p)
  
  # Network information
  Ne <- rowSums(A)
  Ne_treated <- as.vector(A %*% W)
  G = rep(1, N)
  G[Ne_treated == 0] <- 0

  

  # Remove isolates
  W <- W[Ne > 0]
  G <- G[Ne > 0]
  M <- as.numeric(M[Ne > 0])
  X <- X[Ne > 0, ]
  p <- p[Ne > 0]
  N <- length(W)
  A <- A[Ne > 0, Ne > 0]

  # Generate Potential Outcomes
  if (het) {
    x1 <- X[,1]
    tau <- rep(0, N)
    tau[x1==0] <- h
    tau[x1==1] <- - h

    # Generate Treatment Effects
    y0 <- rnorm(N, sd = 0.01)
    y1 <- y0 + tau
    # Generate Outcome
    Y <- y0 * (1-W) + y1 * W
  } else {
    tau <- rep(h, N)
    # Generate Treatment Effects
    y0 <- rnorm(N, sd = 0.01)
    y1 <- y0 + tau
    # Generate Outcome
    Y <- y0 * (1-W) +  y1 * W
  }

  dataset <- list(X = X,
                  Y = Y,
                  W = W,
                  A = A,
                  G = G,
                  M = M,
                  p = p)
  return(dataset)
}
