#' @title
#' Individual Marginal Exposure Probabilities
#'
#' @description
#' Computes the individual marginal probability to be exposed to a given
#' level of the joint intervention.
#' Assumes independent Bernoulli treatment assignment, potentially with
#' heterogeneous unit-level assignment probabilities.
#' Neighborhood exposure is defined as having at least one treated neighbor.
#'
#' @param  i Unit ID
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne_list List of N elements, where each element contains
#' the IDs of the direct neighbors of the corresponding unit
#'
#' @return A numeric value bounded between 0 and 1 measuring
#' the individual marginal probability to be exposed
#' to the level (w,g) of the joint intervention

pi <- function(i, w, g, p, Ne_list) {
  
  vapply(i, function(ii) {
    
    neighbors <- Ne_list[[ii]]
    
    prob_g0 <- prod(1 - p[neighbors])
    prob_g1 <- prob_any_treated(neighbors, p)
    
    prob_w <- (p[ii] ^ w) * ((1 - p[ii]) ^ (1 - w))
    
    prob_w *
      (prob_g1 ^ g) *
      (prob_g0 ^ (1 - g))
  }, numeric(1))
}

# Numerically stable probability that at least one neighbor is treated.
prob_any_treated <- function(neighbors, p) {
  -expm1(sum(log1p(-p[neighbors])))
}
