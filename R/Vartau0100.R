#' @title
#' Estimated Variance for the effect 0100
#'
#' @description
#' Computes the estimated variance of the effect 0100. This effect
#'  compares individuals with the joint treatment exposure (w=0,g=1)
#' with those with (w=0,g=0)
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
Vartau0100 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(0,1)
  vary01 <- Vary(w = 0, g = 1, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(0,0)
  vary00 <- Vary(w = 0, g = 0, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Covariance of the Average Potential Outcomes Y(w,g)=Y(0,1) and Y(w,g)=Y(0,0) 
  covy01y00 <- Covy(w1 = 0, g1 = 1, w2 = 0, g2 = 0,
                    N = N, W = W, G = G, p = p,
                    Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  if (any(is.na(c(vary01, vary00)))) {
    var0100 = NA
  } else {
    var0100 = abs((1 / (N ^ (2))) *
                    (vary01 + vary00 - 2 * covy01y00))
  }
  return(var0100)
}
