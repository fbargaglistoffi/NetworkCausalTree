#' @title
#' Simulate a Clustered Network Environment

#' @description
#' Generates a simulated clustered network environment
#'
#' @param  m Number of clusters
#' @param  N Number of units
#' @param  method_networks method to generate the m networks: "ergm" (Exponential Random Graph Models) ,
#'  "er" (Erdos Renyi) ,"sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er If method "er", probability of the ER model
#' @param  var_homophily_ergm  Variable to account for homophily
#' @param  coef_ergm If method "ergm", coefficients of the ERGM model
#'
#' @return: An adjacency matrix which describes a clustered network environment
#'
#'
generate_clustered_networks = function(m, N, method_networks, param_er,
                              var_homophily_ergm, coef_ergm){

  # Initialize
  comba <- matrix(0, N, N)

  # Compute cluster size
  cluster_size <- N / m

  # Generate m networks
  for (i in 1:m) {

    # ERGM networks
    if (method_networks=="ergm") {
      test.net <- network(cluster_size, directed = FALSE, density = 0)
      test.net%v%"x1" = var_homophily_ergm[(cluster_size*i-(cluster_size-1))]
      g <- simulate(test.net ~ nodematch("x1") + edges,
                        coef = coef_ergm)
      comba[(cluster_size * i - (cluster_size - 1)) : (cluster_size * i),
            (cluster_size * i - (cluster_size - 1)) : (cluster_size * i)] <- as.matrix(g)
    }


    if(method_networks=="er" | method_networks=="sf"){

      # Erdos Renyi networks
      if (method_networks=="er") {
        g=erdos.renyi.game(cluster_size, p.or.m=param_er, type = "gnp")}

      # Barabasi-Albert networks
      if (method_networks=="sf") {
        g=barabasi.game(cluster_size)
      }


      adj<-as.matrix(get.adjacency(g))
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
expand.grid.unique <- function(x, y, include.equals = FALSE){

  x <- unique(x)
  y <- unique(y)
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - include.equals)])
    if (length(z)) cbind(x[i], z, deparse.level=0)
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
shared_neigh = function(i, j, Ne_list){

  shared_neighbors <- length(intersect(Ne_list[[i]], Ne_list[[j]]))

  return(shared_neighbors)
}
