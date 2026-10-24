#' @title
#' Estimated Effect 1000
#'
#' @description
#' Computes the estimate of the treatment effect with no spillover exposure. 
#' This effect compares individuals with the joint treatment exposure (w = 1, g = 0)
#' with those with (w = 0, g = 0)
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

EffTau1000 <- function(N, W, G, Y, p, Ne_list, unit_ids = seq_along(W)) {
  
  idx10 <- which(W == 1 & G == 0)
  idx00 <- which(W == 0 & G == 0)
  
  tau1000 <- 1 / N * (
    sum(Y[idx10] / pi(unit_ids[idx10], 1, 0, p, Ne_list)) -
      sum(Y[idx00] / pi(unit_ids[idx00], 0, 0, p, Ne_list))
  )
  
  return(tau1000)
}
