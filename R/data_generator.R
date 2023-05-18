#' @title
#' Synthetic data generator
#'
#' @description
#' Generates Network Causal Tree synthetic data.
#'
#' @param  N Sample size (default: 2000).
#' @param  m Number of clusters (default: 40).
#' @param  K Number of binary regressors (default: 5).
#' @param  p  N x 1 vector, Probability to be assigned to the active individual
#' intervention (default: rep(0.2,2000))
#' @param het TRUE if the treatment effects 1000 and 1101 are heterogeneous with
#' respect to the first regressor (+taui with X1=0, -taui with X0=0), FALSE if
#' constant (+taui) (default: TRUE).
#' @param taui Absolute value of the treatment effects 1000 and 1101
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
#' - NxN adjacency matrix (`adj_matrix`),
#' - Nx1 individual intervention vector (`W`),
#' - Nx1 neighborhood intervention vector (`G`),
#' - Nx1 outcome vector (`Y`),
#' - Nx1 group membership vector (`M`),
#' - Nx1 degree vector (`Ne`),
#' - Nx1 probability to be assigned to the active individual intervention vector
#' (`p`),
#' - NxK covariates matrix (`X`).
#'
#' @export

data_generator = function(N = 2000,
                          m = 40,
                          K = 5,
                          p = rep(0.2,2000),
                          het = TRUE,
                          taui = 2,
                          method_networks = "er",
                          param_er = 0.1,
                          coef_ergm = NULL,
                          var_homophily_ergm = NULL){

  # Generate Covariates
  X <- NULL
  for(k in 1:K){
    x <- rbinom(N,1,0.5)
    X <- cbind(X,x)
    colnames(X)[k] <- paste0(colnames(X)[k],k)
  }

  # Generate m networks
  gsize <- N/m
  adj_matrix <- genmultnet(N=N,
                           m=m,
                           method_networks=method_networks,
                           param_er = param_er,
                           coef_ergm = coef_ergm,
                           var_homophily_ergm = var_homophily_ergm)

  # Group Indicator
  M <- c(rep(1:m,gsize))
  M <- sort(M)
  levels(M) <- c(1:m)
  net <- graph_from_adjacency_matrix(adj_matrix)

  # Randomly assign unit to treatment arms
  treat <- c()
  for(i in 1:N){
    treat[i] <- rbinom(1,1,prob=p[i])
  }

  #Compute number of treated neigh and consequently Gi
  num_tr_neigh <- as.vector(adj_matrix %*% treat)
  neightreat <- rep(1,N) #Gi
  neightreat[which(num_tr_neigh==0)] <- 0

  # Pass to the standard notation
  neigh <- rowSums(adj_matrix)
  w <- treat[neigh>0]
  g <- neightreat[neigh>0]
  M <- as.numeric(M[neigh>0])
  X <- X[neigh>0,]
  p <- p[neigh>0]
  N <- length(w)
  adj_matrix <- adj_matrix[neigh>0,neigh>0]
  neigh_red <- neigh[neigh>0]

  if (het){
    x1 <- X[,1]
    tau <- rep(0, N)
    tau[x1==0] <- taui
    tau[x1==1] <- - taui

    ## Generate Treatment Effects
    y0 <- rnorm(N,sd=0.01)
    y1 <- y0 + tau
    ## Generate Outcome
    y <- y0*(1-w) + y1*w
  } else {
    tau <- rep(taui, N)
    ## Generate Treatment Effects
    y0 <- rnorm(N,sd=0.01)
    y1 <- y0 + tau
    ## Generate Outcome
    y <- y0*(1-w) +  y1*w
  }

  dataset <- list(adj_matrix = adj_matrix,
                  W = w,
                  G = g,
                  Y = y,
                  M = M,
                  Ne = neigh_red,
                  p = p,
                  X = X)
  return(dataset)
}
