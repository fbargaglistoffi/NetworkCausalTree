#' @title
#' Simulate a Clustered Network Environment

#' @description
#' Generates a simulated clustered network environment
#'
#' @param  k Number of clusters
#' @param  N Number of units
#' @param X N x M Observed Covariates Matrix.
#' @param  method_networks method to generate the k networks: "ergm" (Exponential Random Graph Models) ,
#'  "er" (Erdos Renyi) ,"sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er If method "er", probability of the ER model
#' @param  var_homophily_ergm  Variable to account for homophily
#' @param  coef_ergm If method "ergm", coefficients of the ERGM model
#'
#' @return: An adjacency matrix which describes a clustered network environment
#'
#' @importFrom intergraph asNetwork asIgraph
#'
generate_clustered_networks = function(k,
                                       N,
                                       X = NULL,
                                       method_networks,
                                       param_er,
                                       var_homophily_ergm,
                                       coef_ergm
                                       ){

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
        g = igraph::sample_gnp(cluster_size, p = param_er, directed = FALSE)
      }

      # Barabasi-Albert networks
      if (method_networks == "sf") {
        g = igraph::sample_pa(cluster_size, directed = FALSE)
      }


      adj <- as.matrix(igraph::as_adjacency_matrix(g))
      comba[(cluster_size * i - (cluster_size - 1)) : (cluster_size * i),
            (cluster_size * i - (cluster_size - 1)) : (cluster_size * i)] <- adj

    }

  }
  return(comba)
}
