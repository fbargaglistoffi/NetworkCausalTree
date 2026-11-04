#' Value of the Objective Function (OF)
#'
#' @description
#' Computes the measure of the Objective Function
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param total_variance Sum of variance contributions from other terminal
#' leaves, excluding the current node, for the penalized objective
#' @param nleafs - to be included if method = "penalized" - number of leafs
#' @param unit_ids Original network indices of the observations in W, G, and Y
#' @param N_tr Total number of observations in the discovery sample
#' @param N_est Total number of observations in the estimation sample
#'
#' @return A numeric value corresponding to the computed  Objective Function
#'
compute_OF_Value <- function(method, alpha, beta, gamma, delta,
                             N, W, G, Y, p, Ne, Ne_list, population_effects,
                             total_variance, nleafs,
                             unit_ids = seq_along(W),
                             N_tr = N, N_est = NULL) {

  # Initialize
  inof=NULL

  # Composite criterion
  if (method == "composite") {
    
    if (!is.numeric(population_effects) || length(population_effects) != 4) {
      stop("population_effects must be a numeric vector of length 4")
    }
    if (any(!is.finite(population_effects))) {
      stop("population_effects contains non-finite values: cannot normalize composite OF")
    }
    if (any(population_effects == 0, na.rm = TRUE)) {
      stop("population_effects contains zero: cannot normalize composite OF")
    }

    inof <- alpha * ((EffTau1000(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[1]) ^ 2) +
            beta * ((EffTau1101(N = N, W = W, G = G, Y = Y,
                                p = p, Ne_list = Ne_list,
                                unit_ids = unit_ids) ^ 2) /
                      (population_effects[2]) ^ 2) +
            gamma * ((EffTau1110(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[3]) ^ 2) +
            delta * ((EffTau0100(N = N, W = W, G = G, Y = Y,
                                 p = p, Ne_list = Ne_list,
                                 unit_ids = unit_ids) ^ 2) /
                       (population_effects[4]) ^ 2)
  } else if (method == "penalized") {
    
    if (is.null(N_est) || N_tr <= 0 || N_est <= 0) {
      stop("Penalized objective requires positive N_tr and N_est.")
    }
    
    penalty_coeff <- 1 / N_tr + 1 / N_est

    reward <- compute_OF_Value(
      method = "singular",
      alpha = alpha, beta = beta,
      gamma = gamma, delta = delta,
      N = N, W = W, G = G, Y = Y,
      p = p, Ne = Ne, Ne_list = Ne_list,
      population_effects = population_effects,
      total_variance = NULL, nleafs = nleafs,
      unit_ids = unit_ids
    )
    
    if (alpha > 0) {
      variance <- Vartau1000(
        N = N, W = W, G = G, Y = Y,
        p = p, Ne = Ne, Ne_list = Ne_list,
        unit_ids = unit_ids
      )
    } else if (beta > 0) {
      variance <- Vartau1101(
        N = N, W = W, G = G, Y = Y,
        p = p, Ne = Ne, Ne_list = Ne_list,
        unit_ids = unit_ids
      )
    } else if (gamma > 0) {
      variance <- Vartau1110(
        N = N, W = W, G = G, Y = Y,
        p = p, Ne = Ne, Ne_list = Ne_list,
        unit_ids = unit_ids
      )
    } else if (delta > 0) {
      variance <- Vartau0100(
        N = N, W = W, G = G, Y = Y,
        p = p, Ne = Ne, Ne_list = Ne_list,
        unit_ids = unit_ids
      )
    } else {
      stop("Penalized objective requires a selected treatment effect.")
    }
    
    if (is.null(total_variance)) {
      total_variance <- 0
    }
    
    inof <- reward - penalty_coeff * (total_variance + variance)
    
  } else if (method == "singular") {

   inof <- alpha * (EffTau1000(N = N, W = W, G = G, Y = Y,
                               p = p, Ne_list = Ne_list,
                               unit_ids = unit_ids) ^ 2) +
           beta * (EffTau1101(N = N, W = W, G = G, Y = Y,
                              p = p, Ne_list = Ne_list,
                              unit_ids = unit_ids) ^ 2) +
           gamma * (EffTau1110(N = N, W = W, G = G, Y = Y,
                               p = p, Ne_list = Ne_list,
                               unit_ids = unit_ids) ^ 2) +
           delta * (EffTau0100(N = N, W = W, G = G, Y = Y,
                               p = p, Ne_list = Ne_list,
                               unit_ids = unit_ids) ^ 2)
  }
  
  if (is.null(inof)){
    stop(paste("Unknown method:", method))  
  }

  return(inof)
}
