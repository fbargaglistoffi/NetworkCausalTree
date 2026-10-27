#' @title
#' Estimated Covariance of Average Potential Outcomes
#'
#' @description
#' Computes the estimated covariance of the average potential outcomes
#' related to two given levels of the joint intervention
#' - denoted with (w1,g1) and (w2,g2), respectively -.
#'
#' @param  N Sample size
#' @param  w1 Individual Treatment level - 1
#' @param  g1 Neighborhood Treatment level - 1
#' @param  w2 Individual Treatment level - 2
#' @param  g2 Neighborhood Treatment level -2
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param unit_ids Original network indices of the units included
#' in the covariance calculation
#'
#' @return A numeric value corresponding to the estimated covariance of the
#' average potential outcome under the level w1,g1 and w2,g2 of the joint intervention

Covy = function(N, w1, g1, w2, g2,
                W, G, Y, p,
                Ne, Ne_list, unit_ids = seq_along(W)) {
  
  if (length(W) > 0) {
    
    idx1 <- which(W == w1 & G == g1)
    idx2 <- which(W == w2 & G == g2)
    
    pair_list <- list()
    
    for (i in idx1) {
      id_i <- unit_ids[i]
      
      direct <- Ne_list[[id_i]]
      indirect <- unique(unlist(Ne_list[direct], use.names = FALSE))
      
      candidates <- unique(c(direct, indirect))
      j_candidates <- idx2[unit_ids[idx2] %in% candidates]
      
      for (j in j_candidates) {
        pair_list[[length(pair_list) + 1L]] <- c(i, j)
      }
    }
    
    pairs <- if (length(pair_list) > 0) {
      as.data.frame(t(do.call(cbind, pair_list)))
    } else {
      data.frame(V1 = integer(0), V2 = integer(0))
    }
  
  variab = numeric(nrow(pairs))
    
  for (k in seq_len(nrow(pairs))) {
      
      i <- pairs[k, 1]
      j <- pairs[k, 2]
      
      id_i <- unit_ids[i]
      id_j <- unit_ids[j]
      
      if (shared_neigh(id_i, id_j, Ne_list = Ne_list) > 0 ||
          id_j %in% Ne_list[[id_i]] ||
          id_i %in% Ne_list[[id_j]]) {
        
        pij_val <- pij(id_i, id_j, w1, g1, w2, g2, Ne_list = Ne_list, p = p)
        if (pij_val == 0) {
          next
        }
        pi_i <- pi(id_i, w1, g1, p, Ne_list)
        pi_j <- pi(id_j, w2, g2, p, Ne_list)
        
        variab[k] <- ((pij_val - pi_i * pi_j) / pij_val) *
          (Y[i] / pi_i) * (Y[j] / pi_j)
      }
    }
  
    zero_correction <- 0
    
    for (i in seq_along(W)) {
      id_i <- unit_ids[i]
      
      direct <- Ne_list[[id_i]]
      indirect <- unique(unlist(Ne_list[direct], use.names = FALSE))
      candidates <- unique(c(direct, indirect))
      
      # Include the diagonal pair and dependent units in the leaf
      j_candidates <- unique(c(
        i,
        which(unit_ids %in% candidates)
      ))
      
      for (j in j_candidates) {
        if (i == j) {
          if (W[i] == w1 && G[i] == g1) {
            zero_correction <- zero_correction +
              Y[i]^2 / (2 * pi(unit_ids[i], w1, g1, p, Ne_list))
          }
          
          if (W[j] == w2 && G[j] == g2) {
            zero_correction <- zero_correction +
              Y[j]^2 / (2 * pi(unit_ids[j], w2, g2, p, Ne_list))
          }
          
          next
        }
        
        id_i <- unit_ids[i]
        id_j <- unit_ids[j]
        
        if (shared_neigh(id_i, id_j, Ne_list = Ne_list) == 0 &&
            !(id_j %in% Ne_list[[id_i]]) &&
            !(id_i %in% Ne_list[[id_j]])) {
          next
        }
        
        joint_prob <- pij(
          unit_ids[i], unit_ids[j],
          w1, g1, w2, g2,
          Ne_list = Ne_list, p = p
        )
        
        if (joint_prob == 0) {
          if (W[i] == w1 && G[i] == g1) {
            zero_correction <- zero_correction +
              Y[i]^2 / (2 * pi(unit_ids[i], w1, g1, p, Ne_list))
          }
          
          if (W[j] == w2 && G[j] == g2) {
            zero_correction <- zero_correction +
              Y[j]^2 / (2 * pi(unit_ids[j], w2, g2, p, Ne_list))
          }
        }
      }
    }
    
    covy <- sum(variab) - zero_correction
    
    } else {covy=NA}
  
  return(covy)
}
