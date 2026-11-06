#' @title
#' Split Objective Function
#'
#' @description
#' Splits the sample where the  Objective Function is maximized
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param X N x M matrix, Observed Covariate Matrix
#' @param p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param total_variance - to be included if method = "penalized" - whole variance
#' @param nleafs - to be included if method = "penalized" - number of leafs
#' @param unit_ids Original network indices of the observations in W, G, and Y
#' @param N_tr Total number of observations in the discovery sample
#' @param N_est Total number of observations in the estimation sample
#'
#' @return A numeric vector made up by three elements: the first one identifies the
#' value of the  Objective Function is max , the second one reports the value of the variable that maximizes the  Objective Function,
#' the third one reports the corresponding variable
#'
compute_OF_Split <- function(method, alpha, beta, gamma, delta,
                             W, G, Y, X, p, Ne, Ne_list,
                             population_effects, total_variance, nleafs,
                             unit_ids = seq_along(W),
                             N_tr = length(W), N_est = NULL) {
  
  if (method == "penalized" &&
      (is.null(N_est) || N_tr <= 0 || N_est <= 0)) {
    stop("Penalized splitting requires positive N_tr and N_est.")
  }
  
  # Initialize
  of<- c()
  name <-c()
  values<-c()
  
  if (method == "penalized") {
    parent_reward <- compute_OF_Value(
      method = "singular",
      alpha = alpha, beta = beta,
      gamma = gamma, delta = delta,
      N = length(W),
      W = W, G = G, Y = Y,
      p = p, Ne = Ne, Ne_list = Ne_list,
      unit_ids = unit_ids,
      population_effects = population_effects,
      total_variance = NULL,
      nleafs = nleafs
    )
  }
  
  # Loop over all the predictors, identify all the unique splits for each variable,
  # and compute the corresponding value of OF
  
  for (j in 1:dim(X)[2]) {
    
    x = X[,j]
    
    # sorted unique values
    ux <- sort(unique(x[!is.na(x)]))

    if (length(ux) <= 1) next
    splits <- ux[-1]
    
    ofx <- c()
    
    # Loop over all the possible splits
    
    for (i in seq_along(splits)) {
      
      sp <- splits[i]

      idx_left  <- !is.na(x) & x <  sp
      idx_right <- !is.na(x) & x >= sp
      
      if (!all(xor(idx_left, idx_right))) {
        ofx[i] <- NA
        next
      }
      
      left_ids <- unit_ids[idx_left]
      right_ids <- unit_ids[idx_right]
      
      left_table  <- table(factor(W[idx_left],  levels = c(0,1)), factor(G[idx_left],  levels = c(0,1)))
      right_table <- table(factor(W[idx_right], levels = c(0,1)), factor(G[idx_right], levels = c(0,1)))
      
      if (any(left_table < 1) || any(right_table < 1)) {
        ofx[i] <- NA
        next
      }
      
      n_left <- sum(idx_left)
      n_right <- sum(idx_right)
      n_total <- n_left + n_right
      
      if (method == "penalized") {
        
        left_reward <- compute_OF_Value(
          method = "singular",
          alpha = alpha, beta = beta, gamma = gamma, delta = delta,
          N = n_left, W = W[idx_left], G = G[idx_left],
          Y = Y[idx_left], p = p, Ne = Ne, Ne_list = Ne_list,
          unit_ids = left_ids,
          population_effects = population_effects,
          total_variance = NULL, nleafs = nleafs
        )
        
        right_reward <- compute_OF_Value(
          method = "singular",
          alpha = alpha, beta = beta, gamma = gamma, delta = delta,
          N = n_right, W = W[idx_right], G = G[idx_right],
          Y = Y[idx_right], p = p, Ne = Ne, Ne_list = Ne_list,
          unit_ids = right_ids,
          population_effects = population_effects,
          total_variance = NULL, nleafs = nleafs
        )
        
        reward_change <- (n_left / n_total) * left_reward +
          (n_right / n_total) * right_reward - parent_reward
        
        reward_change <- (n_total / N_tr) * reward_change
        
      }
      
      if (method == "penalized") {
        weighted_variance <- function(ids) {
          n <- length(ids)
          
          args <- list(
            N = n,
            W = W[ids],
            G = G[ids],
            Y = Y[ids],
            p = p,
            Ne = Ne,
            Ne_list = Ne_list,
            unit_ids = unit_ids[ids]
          )
          
          if (alpha > 0) {
            return(do.call(Vartau1000, args))
          }
          if (beta > 0) {
            return(do.call(Vartau1101, args))
          }
          if (gamma > 0) {
            return(do.call(Vartau1110, args))
          }
          if (delta > 0) {
            return(do.call(Vartau0100, args))
          }
          
          stop("No treatment effect selected.")
        }
        
        variance_change <- weighted_variance(which(idx_left)) +
          weighted_variance(which(idx_right)) -
          weighted_variance(seq_along(W))
        penalty_coeff <- 1 / N_tr + 1 / N_est
        
        ofx[i] <- reward_change - penalty_coeff * variance_change
        
      } else {
      ofx[i] <- (n_left / n_total) *
        compute_OF_Value(method, alpha, beta, gamma, delta,
                         N = n_left, W = W[idx_left], G = G[idx_left],
                         Y = Y[idx_left],
                         Ne = Ne, p = p, Ne_list = Ne_list,
                         unit_ids = left_ids,
                         population_effects = population_effects,
                         nleafs = nleafs, total_variance = total_variance,
                         N_tr = N_tr,
                         N_est = N_est) +
        (n_right / n_total) *
        compute_OF_Value(method, alpha, beta, gamma, delta,
                         N = n_right, W = W[idx_right], G = G[idx_right],
                         Y = Y[idx_right],
                         Ne = Ne, p = p, Ne_list = Ne_list,
                         unit_ids = right_ids,
                         population_effects = population_effects,
                         nleafs = nleafs, total_variance = total_variance,
                         N_tr = N_tr,
                         N_est = N_est)
      }
    }
      
    
    valid_idx <- is.finite(ofx)
    if (sum(valid_idx) == 0) next
    
    ofx <- ofx[valid_idx]
    namex  <- rep(colnames(X)[j], sum(valid_idx))
    valuesx <- splits[valid_idx]
    
    # Append all the computed values of the OF, all the variables that have defined the split
    # and their corresponding value.
    
    of = c(of, ofx)
    name = c(name, namex)
    values = c(values, valuesx)
  }
  
  # Identify the variable and the exact splits which maximizes OF and the corresponding OF value
  
  if (all(is.na(of))) {
    ofres <- NA
    splitres <- NA
    varres <- NA
  } else {
    ofres <- max(na.omit(of))
    splitres <- values[which.max(of)]
    varres <- name[which.max(of)]
  }
  
  return(c(of = ofres , split = splitres , var=varres))
  
}
