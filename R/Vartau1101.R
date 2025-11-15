#' @title
#' Estimated Variance for the effect 1101
#'
#' @description
#' Computes the estimated variance of the effect 1000. This effect
#'  compares individuals with the joint treatment exposure (w=1,g=1)
#' with those with (w=0,g=1)
#'
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.
#'
Vartau1101 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(1,1)  
  vary11 <- Vary(w = 1, g = 1, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(0,1)
  vary01 <- Vary(w = 0, g = 1, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Covariance of the Average Potential Outcomes Y(w,g)=Y(1,1) and Y(w,g)=Y(0,1) 
  covy11y01 <- Covy(w1 = 1, g1 = 1,w2 = 0,g2 = 1,
                    N = N,W = W,G = G,
                    p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  if (any(is.na(c(vary11, vary01)))) {
    var1101 = NA
  } else {
    var1101 = abs((1 / (N ^2)) * 
                    (vary11 + vary01 - 2 * covy11y01))
  }
  return(var1101)
}