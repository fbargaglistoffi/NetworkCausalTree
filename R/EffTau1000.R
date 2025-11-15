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
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'

EffTau1000 = function(N, W, G, Y, p, Ne){
tau1000 = 1 / N * (sum(Y[W==1 & G==0] / pi(which(W==1 & G==0), 1, 0, p, Ne)) -
                   sum(Y[W==0 & G==0] / pi(which(W==0 & G==0), 0, 0, p, Ne)))
return(tau1000)
}