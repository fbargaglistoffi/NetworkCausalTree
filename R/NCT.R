#' @title
#' Network Causal Tree

#' @description
#' Returns a Network Causal Tree, with the corresponding estimates
#'
#' @param X N x K Observed Covariates Matrix.
#' @param Y N x 1 Observed Outcome vector.
#' @param W N x 1 Individual Treatment vector.
#' @param effect_weights Treatment Effect weights vector (4 elements):
#' - alpha: weight associated to the treatment effect 1000 (effect of the
#' individual treatment, with the neighborhood treatment set at 0),
#' - beta: weight associated to the treatment effect 1101 (effect of the
#' individual treatment, with the neighborhood treatment set at 1),
#' - gamma: weight associated to the spillover effect 1110 (effect of the
#' neighborhood treatment, with the individual treatment set at 1),
#' - delta: weight associated to the spillover effect 0100 (effect of the
#' neighborhood treatment, with the individual treatment set at 0).
#' @param A N x N Adjacency matrix.
#' @param M N x 1 Cluster Membership vector.
#' @param p  N x 1 Probability to be assigned to the active individual
#' intervention vector.
#' @param ratio_disc Ratio of clusters to be assigned to the discovery set only.
#' @param minsize Minimum number of observaztions for each level of the joint
#' intervention to be required in the leafs
#' @param depth Depth of the tree.
#' @param method Method to compute the objective function: "singular" for NCT
#' targeted to one single effect; "composite" for NCT targeted to multiple
#' effects; "penalized" for a OF computed while considering a single effect only
#' and including a penalization term related to the variance within the leafs.
#' @param output Desired output of the analysis. if output = "detection" only
#' point estimates are computed, if output = "estimation" both estimated effects
#' and variances are computed
#'
#' @return A data.frame describing the obtained Network Causal Trees.
#' Each row represents a partition (of a specific tree) with 10/14 entries.
#' Columns summary:
#' - `OF`: value of the OF in the corresponding partition,
#' - `NOBS_TR`: number of training observations in the partition,
#' - `FILTER`: values of the covariates `X` that identify the partition,
#' - `NOBS_EST`: number of estimation observations in the partition,
#' - `EFF1000_EST`: estimated 1000 effects in the partitions,
#' - `EFF1101_EST`: estimated 1101 effects in the partitions,
#' - `EFF1110_EST`: estimated 1110 effects in the partitions,
#' - `EFF0100_EST`: estimated 0100 effects in the partitions.
#' Additional columns summary (only if output = "Estimation"):
#' - `SETAU1000_EST`: estimated std. error of the 1000 effect in the partition,
#' - `SETAU1101_EST`: estimated std. error of the 1101 effect in the partition,
#' - `SETAU1110_EST`: estimated std. error of the 1110 effect in the partition,
#' - `SETAU0100_EST`: estimated std. error of the 0100 effect in the partition.
#'
#' @import stringi
#' @import statnet
#' @import intergraph
#' @import ergm
#' @import plyr
#' @import stats
#'
#' @export



NetworkCausalTree <- function(X, Y, W,
                               effect_weights = c(1,0,0,0),
                               A = NULL,
                               M = NULL,
                               p = NULL,
                               ratio_disc,
                               depth = 3,
                               minsize = 10,
                               method = "singular",
                               output = "estimation"){

  # compute sample size and number of clusters
  N <- length(W)
  m <- length(unique(M))

  # get effects - specific input weights
  alpha <- effect_weights[1]
  beta <- effect_weights[2]
  gamma <- effect_weights[3]
  delta <- effect_weights[4]

  # check the validity of input parameters
  if (alpha + beta + gamma + delta != 1) {
    stop('Effect weights (effect_weights) must sum up to one.')
  }

  if ((length(which(effect_weights > 0)) > 1) & (method == "singular" | method == "penalized")) {
    stop('If method is set to singular or penalized only one effect should have positive weight.')
  }

  if (1 %in% effect_weights & method == "composite") {
    stop('Composite objective function is computed if at least two effects are investigated.')
  }

  if (is.null(A)) {
    stop('You have to specify the Adiacency Matrix A.')
  }

  if (ratio_disc <= 0 | ratio_disc > 1) {
    stop('The ratio of clusters to be assigned to the discovery set must be above 0 and below 1')
  }

  # compute network - specific information (the degree Ne, the number of treated friends Ne_treated,
  # the value of the neighborhood exposure G and the list of direct neighbors Ne_list)

    Ne <- rowSums(A)
    Ne_treated <- as.vector(A %*% W)
    G = rep(1,N)
    G[Ne_treated == 0] <- 0
    Ne_list <- vector(mode = "list", length = N)

  for (i in  1:N) {
    Ne_list[[i]] <- which(A[i,] > 0)
  }

  # compute the estimated effects in the whole population
  population_effects <- compute_population_effects(N = N, 
                                                   W = W,
                                                   G = G,
                                                   Y = Y, 
                                                   p = p, 
                                                   Ne = Ne)

  # sample the clusters to be assigned to the discovery set
  nclusters_disc = round(m * ratio_disc)
  clusters_disc <- sample(1 : m,
                          size = nclusters_disc, 
                          replace = FALSE)

  # generate the tree on the discovery set
  tree <- sprout_nct(method = method,
                     sampled_clusters = clusters_disc,
                     m = m,
                     depth = depth,
                     minsize = minsize,
                     alpha = alpha,
                     beta = beta,
                     gamma = gamma,
                     delta = delta,
                     N = N,
                     W = W,
                     G = G,
                     Y = Y,
                     X = X,
                     M = M,
                     Ne = Ne,
                     p = p,
                     population_effects = population_effects,
                     Ne_list = Ne_list)


  # organize the tree object
  final_partition <- data.frame(OF = tree$OF,
                                FILTER =  c("NA", as.vector(na.omit(unique(tree$FILTER)))),
                                TERMINAL = tree$TERMINAL,
                                NOBS = tree$NOBS,
                                stringsAsFactors = FALSE)


  # compute the estimates on the partitions generated by the tree
  final_partition_with_estimates <- compute_effects_nct(output = output,
                                                        nct_partition = final_partition,
                                                        N = N,
                                                        W = W,
                                                        G = G,
                                                        Y = Y,
                                                        X = X,
                                                        Ne = Ne,
                                                        p = p,
                                                        Ne_list = Ne_list,
                                                        minsize = minsize)

  return(final_partition_with_estimates)

}
