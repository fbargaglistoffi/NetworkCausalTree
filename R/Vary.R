#' @title
#' Estimated Variance of the Average Potential Outcomes
#'
#' @description
#' Computes the estimated variance of the average potential outcome
#' related to a given level of the joint intervention - denoted with (w,g) -.
#'
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  Y N x 1 vector, Observed Outcome
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the
#'  average potential outcome under the level w,g of the joint intervention


Vary = function(w, g, Y, W, G, 
                p, Ne, Ne_list) {

  if (length(which(W == w & G == g)) > 1) {
    
    pairs <- expand.grid.unique(which(W == w & G == g), 
                                which(W == w & G == g),
                                include.equals = FALSE)
    
    second_element = numeric(nrow(pairs))
    
    for (k in 1 : nrow(pairs)) { 
      
      i = pairs[k, 1]
      j = pairs[k, 2]
      
      if (shared_neigh(i = i, j = j, Ne_list = Ne_list) > 0) {
        
        # covariance correction term
        # accounts for correlation between connected individuals
        pij_val <- pij(i, j, w, g, w, g, Ne, Ne_list, p = p)
        pi_i <- pi(i, w, g, p, Ne)
        pi_j <- pi(j, w, g, p, Ne)
        second_element[k] = ((pij_val - pi_i * pi_j) / pij_val) *
          (Y[i] / pi_i) * (Y[j] / pi_j)
      }
      
      }
    
    # total estimate variance = individual variance component + second_element
    vary = sum((1 - pi(which(W == w & G==g), w, g, p, Ne)) *
           (Y[which(W == w & G == g)] / pi(which(W == w & G == g), w, g, p, Ne)) ^ 2) +
           sum(second_element)
    
  } else {vary = NA}

  return(vary)
}