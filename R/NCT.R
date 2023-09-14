#' @title
#' Network Causal Tree

#' @description
#' Run Network Causal Tree algorithm.
#'
#' @param X N x K Observed Covariates Matrix.
#' @param Y N x 1 Observed Outcome vector.
#' @param W N x 1 Individual Treatment vector.
#' @param A N x N Adjacency matrix.
#' @param M N x 1 Cluster Membership vector.
#' @param p  N x 1 Probability to be assigned to the active individual
#' intervention vector (default: NULL).
#' @param effweights Treatment Effect weights vector (4 elements):
#' - alpha: weight associated to the treatment effect 1000 (effect of the
#' individual treatment, with the neighborhood treatment set at 0),
#' - beta: weight associated to the treatment effect 1101 (effect of the
#' individual treatment, with the neighborhood treatment set at 1),
#' - gamma: weight associated to the spillover effect 1110 (effect of the
#' neighborhood treatment, with the individual treatment set at 1),
#' - delta: weight associated to the spillover effect 0100 (effect of the
#' neighborhood treatment, with the individual treatment set at 0).
#' @param ratio_dis ratio  of clusters to be assigned to the discovery set only
#' (default: 0.5).
#' @param minsize Minimum number of observations for each level of the joint
#' intervention to be required in the leaves (default: 10).
#' @param depth Depth of the tree (default: 3).
#' @param method Method to compute the objective function: "singular" for NCT
#' targeted to one single effect; "composite" for NCT targeted to multiple
#' effects; "penalized" for a OF computed while considering a single effect only
#' and including a penalization term related to the variance
#' (default: "singular").
#' @param output Desired output of the analysis. If output = "detection" only
#' point estimates are computed, if output = "estimation" both estimated effects
#' and variances are computed (default: "estimation").
#'
#' @return A data.frame describing the obtained Network Causal Trees.
#' Each row represents a partition (of a specific tree) with 8/12 entries.
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
#' @import ergm
#' @import plyr
#' @import stats
#'
#' @export
#'
NetworkCausalTrees <- function(X, Y, W, A, M,
                               p = NULL,
                               effweights = c(1,0,0,0),
                               ratio_dis = 0.5,
                               depth = 3,
                               minsize = 10,
                               method = "singular",
                               output = "estimation"){

  N <- length(W)
  m <- length(unique(M))

  # get input weights
  alpha <- effweights[1]
  beta <- effweights[2]
  gamma <-effweights[3]
  delta <- effweights[4]

  # check input parameters
  if (alpha+beta+gamma+delta!=1) {
    stop('Effect weights (effweights) have to sum up to one.')
  }

  if ((length(which(effweights>0))>1) & (method=="singular" | method=="penalized")) {
    stop('If method is set to singular only one effect should have positive weight.')
  }

  if (1 %in% effweights & method=="composite") {
    stop('Composite gof is computed if at least two effects are investigated.')
  }

  if (ratio_dis<=0 | ratio_dis>1) {
    stop('The ratio of clusters to be assigned to the discovery set must be above 0 and below or equal 1')
  }

  Ne <- rowSums(A)
  nt <- as.vector(A %*% W)
  G = rep(1,N)
  G[nt==0] <- 0


  Nel <- vector(mode = "list", length = N)
  for (i in  1:N) {
    Nel[[i]] <- which(A[i,]>0)
  }


  Peff <- popeff(N=N, W=W, G=G, Y=Y, p=p, Ne=Ne)
  mdisc <- round(m*ratio_dis)
  sampgroup_train <- sample(1:m, size=mdisc, replace=FALSE)
  tree <- sproutnetctree(method = method,
                          sampgroup = sampgroup_train,
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
                          Peff = Peff,
                          Nel = Nel)[["tree"]]



  Results <- data.frame(OF = tree$OF,
                       FILTER = c("NA",as.vector(na.omit(unique(tree$FILTER)))),
                       TERMINAL = tree$TERMINAL,
                       NOBS = tree$NOBS,
                       stringsAsFactors = FALSE)


  Results_estimates <- alleffect(output = output,
                                 tree_info = Results,
                                 N = N,
                                 W = W,
                                 G = G,
                                 Y = Y,
                                 X = X,
                                 Ne = Ne,
                                 p = p,
                                 Nel = Nel,
                                 minsize = minsize)
  return(Results_estimates)
}
