#' @title
#' Estimated Covariance of Average Potential Outcomes
#'
#' @description
#' Computes the estimated covariance of the average potential outcomes
#' related to two given levels of the joint intervention
#'  - denoted with (w1,g1) and (w2,g2), respectively -.
#'
#' @param  N Sample size
#' @param  w1 Individual Treatment level - 1
#' @param  g1 Neighborhood Treatment level - 1
#' @param  w2 Individual Treatment level - 2
#' @param  g2 Neighborhood Treatment level -2
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated covariance of the
#'  average potential outcome under the level w1,g1 and w2,g2 of the joint intervention

Covy = function(w1, g1, w2, g2, 
                N, Y, W, G, p, 
                Ne, Ne_list) {
  
  varzero = NULL
  variab = c()
  
  if (length(which(W == w1 & G == g1)) > 1 & 
      length(which(W == w2 & G == g2)) > 1) {
    
    pairs <- expand.grid.unique(which(W == w1 & G == g1),
                                which(W == w2 & G == g2),
                                include.equals = FALSE)
    
    for (k in 1:nrow(pairs)) {
      
      i=pairs[k,1]
      j=pairs[k,2]
      
      if (shared_neigh(i, j, Ne_list = Ne_list) > 0) {
        
        variab = c(varzero,
                   sum(1 / pij(i, j, w1, g1, w2, g2, Ne, Ne_list, p = p)*
                         (Y[i] / pi(i, w1, g1, p, Ne)) * 
                         (Y[j] / pi(j, w2, g2, p, Ne))) *
                     (pij(i, j, w1, g1, w2, g2, Ne, Ne_list, p = p) -
                        pi(i, w1, g1, p, Ne) * pi(j, w2, g2, p,Ne)))
        
      }
      
    }
    
    second_element <- sum(((Y[which(W == w1 & G == g1)]) ^ 2) / (2 * 
                                                                   pi(which(W == w1 & G == g1), w1, g1, p, Ne))) +
      sum(((Y[which(W == w2 & G == g2)]) ^ 2) / (2 * 
                                                   pi(which(W == w2 & G == g2), w2, g2, p, Ne))) 
    
    
    covy = sum(variab) - second_element
    
  } else {covy=NA}
}