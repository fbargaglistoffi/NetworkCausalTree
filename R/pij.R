#' @title
#' Joint Exposure Probabilities
#'
#' @description
#' Computes the individual joint probability of a given pair of units (i and j)
#' to be simultaneously exposed to certain levels of the joint intervention 
#' - (wi, gi) and (wj, gj), respectively -
#' Assumes independent Bernoulli treatment assignment, potentially with
#' heterogeneous unit-level assignment probabilities.
#' Neighborhood exposure is defined as having at least one treated neighbor.
#' Joint probabilities account for shared neighbors and direct connections
#' between units.
#'
#' @param  i Unit ID
#' @param  j Unit ID
#' @param  wi Individual Treatment level of unit i
#' @param  gi Neighborhood Treatment level of unit i
#' @param  wj Individual Treatment level of unit j
#' @param  gj Neighborhood Treatment level of unit j
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#'
#' @return A numeric value bounded between 0 and 1 measuring
#' the individual joint probability of units i and j to be simultaneously
#' exposed to the level wi,gi and wj,gj of  the joint intervention, respectively

pij <- function(i, j, wi, gi, wj, gj, Ne_list, p){
  
  if ((j %in% Ne_list[[i]] && wj == 1 && gi == 0) ||
      (i %in% Ne_list[[j]] && wi == 1 && gj == 0)) {
    return(0)
  }
  
  prob_wi <- (p[i] ^ wi) * ((1 - p[i]) ^ (1 - wi))
  prob_wj <- (p[j] ^ wj) * ((1 - p[j]) ^ (1 - wj))
  
  prob_w <- prob_wi * prob_wj
  
  neighbors_i <- Ne_list[[i]][Ne_list[[i]] != j]
  neighbors_j <- Ne_list[[j]][Ne_list[[j]] != i]
  
  if (setequal(neighbors_i, neighbors_j) &&
      !(j %in% Ne_list[[i]]) &&
      !(i %in% Ne_list[[j]]) &&
      gi != gj) {
    return(0)
  }
  
  all_neighbors <- unique(c(neighbors_i, neighbors_j))
  prob_g00 <- prod(1 - p[all_neighbors])
  
  if (j %in% Ne_list[[i]] && wj == 1) {
    prob_g00 <- 0
  }
  
  if (i %in% Ne_list[[j]] && wi == 1) {
    prob_g00 <- 0
  }
  
  prob_gi0 <- prod(1 - p[neighbors_i])
  prob_gj0 <- prod(1 - p[neighbors_j])
  
  if (j %in% Ne_list[[i]] && wj == 1) {
    prob_gi0 <- 0
  }
  
  if (i %in% Ne_list[[j]] && wi == 1) {
    prob_gj0 <- 0
  }
  
  if (gi == 0 && gj == 0) {
    prob_g <- prob_g00
  } else if (gi == 0 && gj == 1) {
    extra_j <- setdiff(neighbors_j, neighbors_i)
    if (i %in% Ne_list[[j]] && wi == 1) {
      prob_g <- prob_gi0
    } else {
      prob_g <- prob_gi0 * prob_any_treated(extra_j, p)
    }
  } else if (gi == 1 && gj == 0) {
    extra_i <- setdiff(neighbors_i, neighbors_j)
    
    if (j %in% Ne_list[[i]] && wj == 1) {
      prob_g <- prob_gj0
    } else {
      prob_g <- prob_gj0 * prob_any_treated(extra_i, p)
    }
  } else {
    common <- intersect(neighbors_i, neighbors_j)
    only_i <- setdiff(neighbors_i, neighbors_j)
    only_j <- setdiff(neighbors_j, neighbors_i)
    
    q_common <- prod(1 - p[common])
    
    i_has_treated_j <- j %in% Ne_list[[i]] && wj == 1
    j_has_treated_i <- i %in% Ne_list[[j]] && wi == 1
    
    if (i_has_treated_j && j_has_treated_i) {
      prob_g <- 1
    } else if (i_has_treated_j) {
      prob_g <- prob_any_treated(neighbors_j, p)
    } else if (j_has_treated_i) {
      prob_g <- prob_any_treated(neighbors_i, p)
    } else {
      prob_g <- prob_any_treated(common, p) +
        q_common *
        prob_any_treated(only_i, p) *
        prob_any_treated(only_j, p)
    }
  }
  
  pij_val <- prob_w * prob_g
  
  return(pij_val)
  
}
