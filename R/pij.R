#' @title
#' Joint Exposure Probabilities
#'
#' @description
#' Computes the individual joint probability of a given pair of units (i and j)
#'  to be simultaneously exposed to certain levels of the joint intervention 
#'  - (wi, gi) and (wj, gj), respectively -
#'
#' @param  i Unit ID
#' @param  j Unit ID
#' @param  wi Individual Treatment level of unit i
#' @param  wj Individual Treatment level of unit j
#' @param  gi Neighborhood Treatment level of unit i
#' @param  gj Neighborhood Treatment level of unit j
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value bounded between 0 and 1 measuring
#' the individual joint probability of units i and j to be simultaneously
#' exposed to the level wi,gi and wj,gj of  the joint intervention, respectively

pij <- function(i, j, wi, wj, gi, gj, Ne, Ne_list, p){
  pij =  ((p[i] ^ wi) * (1 - p[i]) ^ (1 - wi)) *
    ((p[j] ^ wj) * (1 - p[j]) ^ (1 - wj)) *
    ((1 - (1 - p[i]) ^ (Ne[i] - ifelse (Ne[i] >= Ne[j] & Ne[i] > shared_neigh(i, j, Ne_list = Ne_list),
                                        shared_neigh(i, j, Ne_list = Ne_list),0))) ^ gi) *
    ((1- p[i]) ^ (Ne[i] - ifelse (Ne[i] >= Ne[j] & Ne[i] > shared_neigh(i, j, Ne_list = Ne_list),
                                  shared_neigh(i, j, Ne_list = Ne_list), 0))) ^ (1 - gi)*
    ((1 - (1 - p[j]) ^ (Ne[j] - ifelse (Ne[i] < Ne[j], shared_neigh(i, j, Ne_list = Ne_list),0))) ^ gj) *
    ((1 - p[j])^(Ne[j] - ifelse(Ne[i] < Ne[j], shared_neigh(i, j, Ne_list = Ne_list),0)) ^ (1 - gj))
  return(pij)
}
