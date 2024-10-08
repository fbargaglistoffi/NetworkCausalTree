#' @title
#' Simulate a Clustered Network Environment

#' @description
#' Generates a simulated clustered network environment
#'
#' @param  k Number of clusters
#' @param  N Number of units
#' @param X N x M Observed Covariates Matrix.
#' @param  method_networks method to generate the m networks: "ergm" (Exponential Random Graph Models) ,
#'  "er" (Erdos Renyi) ,"sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er If method "er", probability of the ER model
#' @param  var_homophily_ergm  Variable to account for homophily
#' @param  coef_ergm If method "ergm", coefficients of the ERGM model
#'
#' @return: An adjacency matrix which describes a clustered network environment
#'
#' @import intergraph
#'
generate_clustered_networks = function(k,
                                       N,
                                       method_networks,
                                       param_er,
                                       var_homophily_ergm,
                                       coef_ergm,
                                       X = NULL){

  # Initialize
  comba <- matrix(0, N, N)

  # Compute cluster size
  cluster_size <- N / k

  # Generate k networks
  for (i in 1:k) {

    # ERGM networks
    if (method_networks == "ergm") {
      test.net <- make_empty_graph(n = cluster_size,
                                   directed = FALSE)
      test.net <- intergraph::asNetwork(test.net)
      homophily_vector <- X[,var_homophily_ergm]
      test.net%v%"homophily" =  homophily_vector[(cluster_size * i -
                                                 (cluster_size - 1)) :
                                                 (cluster_size * i) ]
      g <- ergm::simulate_formula(test.net ~ nodematch("homophily") + edges,
                                  coef = coef_ergm)
      comba[(cluster_size * i - (cluster_size - 1)) : (cluster_size * i),
            (cluster_size * i - (cluster_size - 1)) : (cluster_size * i)] <- as.matrix(g)
    }


    if (method_networks =="er" | method_networks == "sf") {

      # Erdos Renyi networks
      if (method_networks == "er") {
        g = igraph::erdos.renyi.game(cluster_size,
                                   p.or.m=param_er,
                                   type = "gnp")}

      # Barabasi-Albert networks
      if (method_networks == "sf") {
        g = igraph::barabasi.game(cluster_size)
      }


      adj <- as.matrix(igraph::as_adjacency_matrix(g))
      comba[(cluster_size * i - (cluster_size - 1)) : (cluster_size * i),
            (cluster_size * i - (cluster_size - 1)) : (cluster_size * i)] <- adj

    }

  }
  return(comba)
}


#' @title
#' Generate all combinations of the supplied vectors, without rows characterized
#' by equal elements
#' @description
#' Create a data frame including all combinations of the supplied vectors or factors,
#' while omitting rows characterized by equal elements.
#' @param  x A vector
#' @param  y A vector
#' @param  include.equals Boolean (dafault: FALSE)
#'
#' @return A data frame with all combinations of the elements of the vectors,
#' with no rows characterized by equal elements
#'
#' @import dplyr
#'
expand.grid.unique <- function(x, y, include.equals = FALSE){

  x <- unique(x)
  y <- unique(y)

  g <- function(i){
    z <- dplyr::setdiff(y, x[seq_len(i - include.equals)])
    if (length(z)) cbind(x[i], z, deparse.level = 0)
  }

  do.call(rbind, lapply(seq_along(x), g))

}



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
