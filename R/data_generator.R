#' @title
#' Synthetic data generator
#'
#' @description
#' Generates Network Causal Tree synthetic data.
#'
#' @param  N Sample size (default: 2000).
#' @param  M Number of binary regressors (default: 5).
#' @param  k Number of clusters (default: 40).
#' @param  p  N x 1 vector, Probability to be assigned to the active individual
#' intervention (default: rep(0.2,2000))
#' @param h Absolute value of the treatment effects 1000 and 1101
#' (default: 2).
#' @param het TRUE if the treatment effects 1000 and 1101 are heterogeneous with
#' respect to the first regressor (h with X1=0, -taui with X0=0), FALSE if
#' constant (-h) (default: TRUE).
#' @param  method_networks Method to generate the m networks:
#' "ergm" (Exponential Random Graph Models), "er" (Erdos Renyi), "sf"
#' (Barabasi-Albert model) (default: "er").
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er Probability of the "er" model, if used (default: 0.2).
#' @param  coef_ergm Coefficients of the "ergm" model, if used (default: NULL).
#' @param  var_homophily_ergm Variable to account for homophily in the "ergm"
#' @param remove_isolates Logical; remove isolated nodes? (default TRUE)
#' model (default: NULL).
#'
#' @return A list of synthetic data containing:
#' - NxM covariates matrix (`X`).
#' - Nx1 outcome vector (`Y`),
#' - Nx1 individual intervention vector (`W`),
#' - NxN adjacency matrix (`A`),
#' - Nx1 neighborhood intervention vector (`G`),
#' - Nx1 group membership vector (`K`),
#' - Nx1 probability to be assigned to the active individual intervention vector
#' (`p`),
#'
#' @importFrom igraph graph_from_data_frame V E make_empty_graph layout_as_tree "E<-" "V<-"
#'
#' @export

data_generator = function(N = 2000,
                          M = 5,
                          k = 40,
                          p = rep(0.2,2000),
                          h = 2,
                          het = TRUE,
                          method_networks = "er",
                          param_er = 0.2,
                          coef_ergm = NULL,
                          var_homophily_ergm = NULL,
                          remove_isolates = TRUE){

  # check the validity of input parameters
  if (length(p) != N) {
    stop('The length of vector describing individual probabilities to be assigned to the intervention MUST be equal to N')
  }


  # Generate Covariates
  X <- NULL
  for (m in 1 : M) {
    x <- rbinom(N, 1, 0.5)
    X <- cbind(X, x)
    colnames(X)[m] <- paste0(colnames(X)[m], m)
  }

  # Generate m networks
  A <- generate_clustered_networks(N = N,
                                   k = k,
                                   method_networks = method_networks,
                                   param_er = param_er,
                                   coef_ergm = coef_ergm,
                                   var_homophily_ergm = var_homophily_ergm,
                                   X = X)

  net <- igraph::graph_from_adjacency_matrix(A)

  # Group Indicator
  cluster_size <- N / k
  K <- c(rep(1 : k, cluster_size))
  K <- sort(K)
  levels(K) <- c(1 : k)


  # Randomly assign unit to treatment arms
  W <- rbinom(N, 1, prob = p)

  # Network information
  Ne <- rowSums(A)
  Ne_treated <- as.vector(A %*% W)
  G = rep(1, N)
  G[Ne_treated == 0] <- 0



  if (remove_isolates) {
    W <- W[Ne > 0]
    G <- G[Ne > 0]
    K <- as.numeric(K[Ne > 0])
    X <- X[Ne > 0, ]
    p <- p[Ne > 0]
    N <- length(W)
    A <- A[Ne > 0, Ne > 0]
  }

  # Generate Potential Outcomes
  # effect size depends on first covariate
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
                  K = K,
                  p = p)
  return(dataset)
}
