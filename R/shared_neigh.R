#' @title
#' Number of Shared Neighbors
#'
#' @description
#' Computes the number of shared neighbors between unit i and unit j.
#'
#' @param  i Unit ID
#' @param  j Unit ID
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value representing the number of shared neighbors between
#' unit i and j.
#'
#' @import dplyr
#'
#'
shared_neigh = function(i, j, Ne_list){
  
  shared_neighbors <- length(dplyr::intersect(Ne_list[[i]], Ne_list[[j]]))
  
  return(shared_neighbors)
}
