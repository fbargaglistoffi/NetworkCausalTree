#' @title
#' Simulate a Clustered Network Environment

#' @description
#' Generates a simulated clustered network environment.
#'
#' @param m Number of clusters.
#' @param N Number of units.
#' @param method_networks method to generate the m networks: "ergm" (Exponential
#' Random Graph Models) , "er" (Erdos Renyi), "sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a
#' multiple of m.
#' @param param_er If method "er", probability of the ER model.
#' @param var_homophily_ergm  Variable to account for homophily.
#' @param coef_ergm If method "ergm", coefficients of the ERGM model.
#'
#' @return An adjacency matrix which describes a clustered network environment.
#'
genmultnet <- function(m,N,method_networks,param_er,var_homophily_ergm,coef_ergm){

  comba <- matrix(0,N,N)
  gsize <- N/m

  for (i in 1:m) {
    if (method_networks=="ergm") {
      test.net <- network(gsize, directed = FALSE, density = 0)
      test.net%v%"x1" <- var_homophily_ergm[(gsize*i-(gsize-1))]
      g.sim <- simulate(test.net ~ nodematch("x1") + edges,
                        coef = coef_ergm)
      comba[(gsize*i-(gsize-1)):(gsize*i),(gsize*i-(gsize-1)):(gsize*i)] <- as.matrix(g.sim)
    }

    if (method_networks=="er" | method_networks=="sf") {
      if (method_networks=="er") {
        g <- erdos.renyi.game(gsize, p=param_er, type = "gnp")
      }
      else if (method_networks=="sf") {
        g <- barabasi.game(gsize)
      }
      adj <- as.matrix(get.adjacency(g))
      comba[(gsize*i-(gsize-1)):(gsize*i),(gsize*i-(gsize-1)):(gsize*i)] <- adj
    }
  }
  return(comba)
}


#' @title
#' Generate all combinations of the supplied vectors, without rows characterized
#' by equal elements
#'
#' @description
#' Create a data frame including all combinations of the supplied vectors or factors,
#' while omitting rows characterized by equal elements.
#'
#' @param x A vector.
#' @param y A vector.
#' @param include.equals Boolean.
#'
#' @return A data frame with all combinations of the elements of the vectors,
#' with no rows characterized by equal elements.
#'
expand.grid.unique <- function(x, y, include.equals = FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)])
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
#' @param i Unit ID.
#' @param j Unit ID.
#' @param Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i.
#'
#' @return A numeric value representing the number of shared neighbors between
#' unit i and j.
#'
sharedn <- function(i,j,Nel) {
  return(length(intersect(Nel[[i]],Nel[[j]])))
}
