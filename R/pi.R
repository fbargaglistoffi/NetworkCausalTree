#' @title
#' Individual Marginal Exposure Probabilities
#'
#' @description
#' Computes the individual marginal probability to be exposed to a given
#' level of the joint intervention.
#'
#' @param  i Unit ID
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value bounded between 0 and 1 measuring
#' the individual marginal probability to be exposed
#' to the level (w,g) of the joint intervention

pi <- function(i, w, g, p, Ne){
  pi <- (p[i] ^ w) *
        (1 - p[i]) ^ (1 - w) *
        ((1 - (1 - p[i]) ^ Ne[i]) ^ g) *
        ((1 - p[i]) ^ Ne[i]) ^ (1 - g)
  return(pi)
}
