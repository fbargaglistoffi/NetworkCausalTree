#' @title
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
#' @param p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param total_variance - to be included if method = "penalized" - whole variance
#' @param nleafs - to be included if method = "penalized" - number of leafs
#'
#' @return A numeric value corresponding to the computed  Objective Function
#'
compute_OF_Value = function(method, alpha, beta, gamma, delta,
                   N, W, G, Y, p, Ne, Ne_list, population_effects,
                   total_variance, nleafs){

  # initialize
  inof=NULL

  # composite criterion
  # divide by population effects to normalize effect 
  if (method == "composite") {

    inof <- alpha * (((EffTau1000(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[1]) ^ 2) +
            beta * (((EffTau1101(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                      (population_effects[2]) ^ 2) +
            gamma * (((EffTau1110(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[3]) ^ 2) +
            delta * (((EffTau0100(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[4]) ^ 2)
  }

  # penalized criterion
  if (method == "penalized") {

  inof <-  alpha * (((EffTau1000(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                    2 / nleafs * sum(c(total_variance, Vartau1000(N = N, W = W, Y = Y,G  = G,
                    p = p, Ne = Ne, Ne_list = Ne_list)))) +
           beta *  (((EffTau1101(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                    2 / nleafs *sum(c(total_variance, Vartau1101(N = N, W = W, Y = Y, G = G,
                    p = p, Ne = Ne, Ne_list = Ne_list)))) +
           gamma * (((EffTau1110( N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                   2 / nleafs * sum(c(total_variance,Vartau1110(N = N, W = W, Y = Y, G = G,
                   p = p, Ne = Ne, Ne_list = Ne_list))))+
           delta * (((EffTau0100(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                   2 / nleafs * sum(c(total_variance,Vartau0100(N = N, W = W, Y = Y, G = G,
                   p = p, Ne = Ne, Ne_list = Ne_list)))
    )
  }

  # singular criterion
  if (method == "singular") {

   inof <- alpha * ((EffTau1000( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2) +
           beta * ((EffTau1101( N = N, W = W, G = G,
                                Y = Y, p = p, Ne = Ne)) ^ 2) +
           gamma * ((EffTau1110( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2) +
           delta * ((EffTau0100( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2)
  }

  return(inof)
}
