#' @title
#' Estimated Effect 1101

#' @description
#' Computes the estimate of the treatment effect with spillover exposure. 
#' This effect compares individuals with the joint treatment exposure (w = 1, g = 1)
#' with those with (w = 0, g = 1)
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne_list List of direct neighbors for each unit in the exposure-defining network
#' @param unit_ids Original network indices of the units included in the estimate
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau1101 <- function(N, W, G, Y, p, Ne_list, unit_ids = seq_along(W)) {
  
  idx11 <- which(W == 1 & G == 1)
  idx01 <- which(W == 0 & G == 1)
  
  tau1101 <- 1 / N * (
    sum(Y[idx11] / pi(unit_ids[idx11], 1, 1, p, Ne_list)) -
      sum(Y[idx01] / pi(unit_ids[idx01], 0, 1, p, Ne_list))
  )
  
  return(tau1101)
}
