#' @title
#' Estimated Effects in the whole sample
#'
#' @description
#' Computes the four Estimated Effects in the Whole Population
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric vector made up by four elements, representing the
#' estimated effects in the whole population: the first element is the effect 1000,
#' the second one refers to the effect 1101, the third one refers to 1110 and the
#' fourth refers to 0100
#'
compute_population_effects = function(N, W, G, Y, p, Ne){
  PEffTau1000 <- EffTau1000(N = N, W = W, G = G , Y = Y, p = p, Ne = Ne)
  PEffTau1101 <- EffTau1101(N = N, W = W, G = G , Y = Y, p = p, Ne = Ne)
  PEffTau1110 <- EffTau1110(N = N, W = W, G = G , Y = Y, p = p, Ne = Ne)
  PEffTau0100 <- EffTau0100(N = N, W = W, G = G , Y = Y, p = p, Ne = Ne)
  Peff <- c(PEffTau1000, PEffTau1101, PEffTau1110, PEffTau0100)
  return(Peff)
}