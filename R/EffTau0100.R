#' @title
#' Estimated Effect 0100
#'
#' @description
#' Computes the estimate of the spillover effect in the absence of an individual 
#' exposure to the intervention. 
#' This effect compares individuals with the joint treatment exposure (w = 0, g = 1)
#' with those with (w = 0, g = 0)
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau0100 = function(N, W, G, Y, p, Ne){
  tau0100 = 1 / N * (sum(Y[W==0 & G==1] / pi(which(W==0 & G==1), 0, 1, p, Ne)) -
                       sum(Y[W==0 & G==0] / pi(which(W==0 & G==0), 0, 0, p, Ne)))
  return(tau0100)
}