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
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param X N x N matrix, Observed Covariate Matrix
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
                            N, W, G, Y, X, p, Ne, Ne_list,
                            population_effects, total_variance, nleafs){
  
  #initialize
  of<- c()
  name <-c()
  values<-c()
  
  # loop over all the predictors, identify all the unique splits for each variable,
  # and compute the corresponding value of OF
  
  for (j in 1:dim(X)[2]) {
    
    x = X[,j]
    splits <- sort(unique(x))
    valuesx<-c()
    ofx <- c()
    namesx <- c()
    
    # loop over all the possible splits
    
    for (i in seq_along(splits)) {
      
      sp <- splits[i]
      
      if (all(as.numeric(table(W[x >= sp], G[x >= sp]))>2) &
          all(as.numeric(table(W[x < sp], G[x < sp])) > 2)) {
        
        # Compute the corresponding OF
        ofx[i]<- 1/2*(compute_OF_Value(method = method, alpha = alpha, beta = beta,
                                       gamma = gamma,delta = delta, N = length(which(x < sp)),
                                       W = W[x < sp], G = G[x < sp],Y = Y[x < sp],
                                       Ne = Ne[x < sp], p = p[x < sp],
                                       population_effects = population_effects,
                                       Ne_list = Ne_list[x < sp],  nleafs = nleafs,
                                       total_variance = total_variance) +
                        compute_OF_Value(method = method, alpha = alpha, beta = beta,
                                         gamma = gamma,delta = delta, N = length(which(x >= sp)),
                                         W = W[x >= sp], G = G[x >= sp], Y = Y[x >= sp],
                                         Ne = Ne[x >= sp], p = p[x >= sp],
                                         population_effects = population_effects,
                                         Ne_list = Ne_list[x >= sp], nleafs = nleafs,
                                         total_variance = total_variance))
      }
      else {ofx[i] <- 0}
    }
    
    namex = rep(colnames(X)[j], length(unique(x)))
    valuesx = c(sort(unique(x)))
    
    # append all the computed values of the OF, all the variables that have defined the split
    # and their corresponding value.
    
    of = c(of, ofx)
    name = c(name, namex)
    values = c(values, valuesx)
  }
  
  # identify the variable and the exact splits which naximizes OF and the corresponding OF value
  
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
