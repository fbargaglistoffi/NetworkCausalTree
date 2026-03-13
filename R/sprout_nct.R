#' @title Generation of the Network Causal Tree
#'
#' @description Sprouts the network causal tree, eventually including a fraction of the initial -discovery-
#' sample.
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
#' @param minsize minimum number of observations for each level of the joint intervention
#' to be required in the leafs
#' @param depth depth of the tree
#'
#' @return A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the OF in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' values of the Xs to identify the given partition, and the column TERMINAL reporting
#' whether the node is a splitting node ("SPLIT") or a leaf ("LEAF").
#'
sprout_nct = function(method, alpha, beta, gamma, delta,
                      N, sampled_clusters, W, G, Y, X, K = NULL, p,
                      Ne, Ne_list, population_effects, minsize, depth){
  
  if (is.null(K)) {
    # If K is not provided, treat all units as belonging to a single cluster
    K <- rep(1, N)
    } else if (length(K) != N) {
    stop("K must be a vector of length N.")
  }
  
  # Initialize data frame
  data <- data.frame(idunit = 1:N, W = W, G = G, Y = Y, K = K)
  
  # Attach X
  if (!is.null(X) && NCOL(X) > 0) {
    # Validate that X has N observations
    if (is.null(dim(X))) {
      # X is a vector
      if (length(X) != N) {
        stop("X must have length N (", N, "), but has length ", length(X), ".")
      }
    } else {
      # X is a matrix/data.frame/array
      if (nrow(X) != N) {
        stop("X must have N rows (", N, "), but has ", nrow(X), " rows.")
      }
    }
    
    X <- as.data.frame(X)
    if (!all(grepl("^X\\.", colnames(X)))) {
      colnames(X) <- paste0("X.", seq_len(ncol(X)))
    }
    data <- cbind(data, X)
  }
  
  if (is.null(sampled_clusters)) {
    # By default, treat all clusters as part of the discovery set
    sampled_clusters <- unique(K)
  }
  # Take only those observations in the discovery set
  datasample <- data[which(K %in% sampled_clusters), ]
  datasample <- datasample[order(datasample$idunit), ]
  sampleid <- unique(datasample$idunit)
  
  if (length(sampleid) == 0) {
    return(data.frame(
      NODE = 1,
      OF = 0,
      NOBS = 0,
      FILTER = NA,
      TERMINAL = "LEAF",
      stringsAsFactors = FALSE
    ))
  }
  
  W = as.numeric(datasample$W)
  G = as.numeric(datasample$G)
  Y = as.numeric(datasample$Y)

  covar_cols <- grep("^X\\.", names(datasample))
  if (length(covar_cols) == 0) {
    return(data.frame(
      NODE = 1,
      OF = 0,
      NOBS = nrow(datasample),
      FILTER = NA,
      TERMINAL = "LEAF",
      stringsAsFactors = FALSE
    ))
  }

  X <- as.matrix(datasample[, covar_cols, drop = FALSE])

  # Identify partitions
  tree_info <- identify_partitions_nct(method = method, alpha = alpha, beta = beta,
                                       gamma = gamma, delta = delta,
                                       N = length(sampleid), depth = depth,
                                       minsize = minsize,
                                       W = W, G = G,Y = Y, X = X,
                                       p = p[sampleid], Ne = Ne[sampleid],
                                       Ne_list = Ne_list[sampleid],
                                       population_effects = population_effects)
  
  return(tree_info)
}
