#' @title
#' Estimated Variance of the Average Potential Outcomes
#'
#' @description
#' Computes the estimated variance of the average potential outcome
#' related to a given level of the joint intervention - denoted with (w,g) -.
#'
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  Y N x 1 vector, Observed Outcome
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param unit_ids Original network indices of the units included
#' in the variance calculation
#'
#' @return A numeric value corresponding to the estimated variance of the
#'  average potential outcome under the level w,g of the joint intervention


Vary = function(w, g, Y, W, G,
                p, Ne, Ne_list, unit_ids = seq_along(W)) {

  if (length(W) > 0) {
    
    idx <- which(W == w & G == g)
    
    pair_list <- list()
    
    for (i in idx) {
      id_i <- unit_ids[i]
      
      direct <- Ne_list[[id_i]]
      indirect <- unique(unlist(Ne_list[direct], use.names = FALSE))
      
      candidates <- unique(c(direct, indirect))
      j_candidates <- idx[unit_ids[idx] %in% candidates]
      j_candidates <- j_candidates[j_candidates != i]
      
      for (j in j_candidates) {
        pair_list[[length(pair_list) + 1L]] <- c(i, j)
      }
    }
    
    pairs <- if (length(pair_list) > 0) {
      as.data.frame(t(do.call(cbind, pair_list)))
    } else {
      data.frame(V1 = integer(0), V2 = integer(0))
    }
    
    second_element = numeric(nrow(pairs))
    
    for (k in seq_len(nrow(pairs))) {
      
      i <- pairs[k, 1]
      j <- pairs[k, 2]
      
      id_i <- unit_ids[i]
      id_j <- unit_ids[j]
      
      if (shared_neigh(id_i, id_j, Ne_list = Ne_list) > 0 ||
          id_j %in% Ne_list[[id_i]] ||
          id_i %in% Ne_list[[id_j]]) {
        
        # covariance correction term
        # accounts for correlation between connected individuals
        pij_val <- pij(id_i, id_j, w, g, w, g, Ne_list = Ne_list, p = p)
        if (pij_val == 0) {
          stop("Zero joint exposure probability encountered in Vary")
        }
        pi_i <- pi(id_i, w, g, p, Ne_list)
        pi_j <- pi(id_j, w, g, p, Ne_list)
        second_element[k] = ((pij_val - pi_i * pi_j) / pij_val) *
          (Y[i] / pi_i) * (Y[j] / pi_j)
      }
      
      }
    
    idx <- which(W == w & G == g)
    prob <- pi(unit_ids[idx], w, g, p, Ne_list)
    
    zero_correction <- 0
    
    pair_list <- list()
    
    for (i in seq_along(W)) {
      id_i <- unit_ids[i]
      
      direct <- Ne_list[[id_i]]
      
      indirect <- unique(unlist(Ne_list[direct], use.names = FALSE))
      
      candidates <- unique(c(direct, indirect))
      
      j_candidates <- which(unit_ids %in% candidates)
      j_candidates <- j_candidates[j_candidates > i]
      
      for (j in j_candidates) {
        pair_list[[length(pair_list) + 1L]] <- c(i, j)
      }
    }
    
    all_pairs <- if (length(pair_list) > 0) {
      do.call(cbind, pair_list)
    } else {
      matrix(integer(0), nrow = 2)
    }
    
    for (k in seq_len(ncol(all_pairs))) {
      i <- all_pairs[1, k]
      j <- all_pairs[2, k]
      
      id_i <- unit_ids[i]
      id_j <- unit_ids[j]
      
      if (shared_neigh(id_i, id_j, Ne_list = Ne_list) == 0 &&
          !(id_j %in% Ne_list[[id_i]]) &&
          !(id_i %in% Ne_list[[id_j]])) {
        next
      }
      
      if (pij(unit_ids[i], unit_ids[j], w, g, w, g,
              Ne_list = Ne_list, p = p) == 0) {
        
        if (W[i] == w && G[i] == g) {
          zero_correction <- zero_correction +
            Y[i]^2 / pi(unit_ids[i], w, g, p, Ne_list)
        }
        
        if (W[j] == w && G[j] == g) {
          zero_correction <- zero_correction +
            Y[j]^2 / pi(unit_ids[j], w, g, p, Ne_list)
        }
      }
    }
    
    vary <- sum((1 - prob) * (Y[idx] / prob)^2) +
      sum(second_element) + zero_correction
    
  } else {vary = NA}

  return(vary)
}