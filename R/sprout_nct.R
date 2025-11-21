#' @title
#' Generation of the Network Causal Tree
#'
#' @description
#
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param  N Sample size
#' @param sampled_clusters Clusters assigned to the discovery set
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  K N x 1 vector, Cluster Membership
#' @param  X N x M matrix, Observed Covariates Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize minimum number of observaztions for each level of the joint intervention
#' to be required in the leafs
#' @param depth depth of the tree
#' @param A Adjacency matrix (internal use)
#'
#' @return A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the OF in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' valyes of the Xs to identify the given partition,, the column TERMINAL reports the
#' 'role' of the node - parent or leaf -.
#'
sprout_nct = function(method,
                      sampled_clusters,
                      depth,
                      minsize,
                      alpha, beta, gamma, delta,
                      N, W, G, Y, X, K,
                      Ne, p, population_effects,
                      Ne_list,
                      A = NULL) {
  
  # create all-zero matrix if no adjacency matrix is passed
  if (is.null(A)) A <- matrix(0, nrow = N, ncol = N)
  
  # compute network features
  Ne <- rowSums(A)
  Ne_treated <- as.vector(A %*% W)
  G <- as.numeric(Ne_treated > 0)
  if (length(G) != N) G <- rep(0, N)
  
  W <- as.numeric(W)[seq_len(N)]
  Y <- as.numeric(Y)[seq_len(N)]
  G <- as.numeric(G)[seq_len(N)]
  K <- as.numeric(K)[seq_len(N)]
  
  df <- data.frame(
    idunit = seq_len(N),
    W = W,
    G = G,
    Y = Y,
    K = K
  )
  
  datasample <- df[K %in% sampled_clusters, , drop = FALSE]
  
  # edge case: no data in discovery sample
  if (nrow(datasample) == 0) {
    return(data.frame(
      NODE = 1,
      OF = NA,
      NOBS = 0,
      FILTER = NA,
      TERMINAL = "LEAF"
    ))
  }
  
  # order and extract discovery subset
  datasample <- datasample[order(datasample$idunit), ]
  sampleid <- datasample$idunit
  W <- datasample$W
  G <- datasample$G
  Y <- datasample$Y
  
  # ensure covariate matrix is consistent with the subset used for discovery
  if (is.null(X) || length(X) == 0) {
    X <- matrix(NA, nrow = nrow(datasample), ncol = 0)
  } else {
    X <- as.matrix(X)[sampleid, , drop = FALSE]
    if (ncol(X) == 0) X <- matrix(NA, nrow = nrow(datasample), ncol = 0)
    if (!is.null(X) && ncol(X) > 0) {
      colnames(X) <- paste0("X.", seq_len(ncol(X)))
    }
  }
  
  # degrees and neighbors for sampled individuals
  Ne_s <- Ne[sampleid]
  Ne_list_s <- Ne_list[sampleid]
  
  tree_info <- identify_partitions_nct(
    method = method,
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    N = length(sampleid),
    depth = depth,
    minsize = minsize,
    W = W, G = G, Y = Y, X = X,
    p = p[sampleid],
    Ne = Ne_s,
    Ne_list = Ne_list_s,
    population_effects = population_effects
  )
  
  return(tree_info)
}