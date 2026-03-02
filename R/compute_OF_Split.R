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
#'
#' @return A numeric vector made up by three elements: the first one identifies the
#' value of the  Objective Function is max , the second one reports the value of the variable that maximizes the  Objective Function,
#' the third one reports the corresponding variable
#'
compute_OF_Split = function(method, alpha, beta, gamma, delta,
                            W, G, Y, X, p, Ne, Ne_list,
                            population_effects, total_variance, nleafs){
  
  # Initialize
  of<- c()
  name <-c()
  values<-c()
  
  # Loop over all the predictors, identify all the unique splits for each variable,
  # and compute the corresponding value of OF
  
  for (j in 1:dim(X)[2]) {
    
    x = X[,j]
    
    # sorted unique values
    ux <- sort(unique(x))

    if (length(ux) > 1) {
      splits <- ux[-1]
    } else {
      splits <- ux
    }
    
    ofx <- c()
    
    # Loop over all the possible splits
    
    for (i in seq_along(splits)) {
      
      sp <- splits[i]

      left_table  <- table(W[!is.na(x) & x < sp],  G[!is.na(x) & x < sp])
      right_table <- table(W[!is.na(x) & x >= sp], G[!is.na(x) & x >= sp])
      
      if (any(left_table < 1) || any(right_table < 1)) {
        ofx[i] <- NA
        next
      }
      
      ofx[i] <- 1/2 * (
        compute_OF_Value(method, alpha, beta, gamma, delta,
                         N = sum(!is.na(x) & x <  sp), W = W[!is.na(x) & x <  sp], G = G[!is.na(x) & x <  sp], Y = Y[!is.na(x) & x <  sp],
                         Ne = Ne[!is.na(x) & x <  sp], p = p[!is.na(x) & x <  sp], Ne_list = Ne_list[!is.na(x) & x <  sp],
                         population_effects = population_effects,
                         nleafs = nleafs, total_variance = total_variance) +
          compute_OF_Value(method, alpha, beta, gamma, delta,
                           N = sum(!is.na(x) & x >= sp), W = W[!is.na(x) & x >= sp], G = G[!is.na(x) & x >= sp], Y = Y[!is.na(x) & x >= sp],
                           Ne = Ne[!is.na(x) & x >= sp], p = p[!is.na(x) & x >= sp], Ne_list = Ne_list[!is.na(x) & x >= sp],
                           population_effects = population_effects,
                           nleafs = nleafs, total_variance = total_variance)
      )
    }
    
    valid_idx <- !is.na(ofx)
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
