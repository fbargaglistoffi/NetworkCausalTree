#' @title
#' Simulate a Clustered Network Environment

#' @description
#' Generates a simulated clustered network environment
#'
#' @param k Number of clusters
#' @param N Number of units
#' @param X N x M Observed Covariates Matrix
#' @param method_networks method to generate the k networks: "ergm" (Exponential Random Graph Models) ,
#' "er" (Erdos Renyi) ,"sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a multiple of k
#' @param param_er If method "er", probability of the ER model
#' @param var_homophily_ergm  Variable to account for homophily
#' @param coef_ergm If method "ergm", coefficients of the ERGM model
#'
#' @return An adjacency matrix which describes a clustered network environment
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
    
    idx <- (cluster_size * i - (cluster_size - 1)):(cluster_size * i)
    
    # ERGM networks
    if (method_networks == "ergm") {
      test.net <- make_empty_graph(n = cluster_size,
                                   directed = FALSE)
      test.net <- intergraph::asNetwork(test.net)
      
      # Extracts the column from X that will be used for homophily
      homophily_vector <- X[,var_homophily_ergm]
      
      # Assign homophily attributes to nodes in this cluster
      test.net%v%"homophily" =  homophily_vector[idx]
      
      # Simulate an ERGM network
      g <- ergm::simulate_formula(test.net ~ nodematch("homophily") + edges,
                                  coef = coef_ergm)
      
      # Places the generated network into the correct block of the matrix
      comba[idx, idx] <- as.matrix(g)
    }


    else if (method_networks =="er" | method_networks == "sf") {

      # Erdos Renyi networks
      if (method_networks == "er") {
        g = igraph::sample_gnp(cluster_size, p = param_er, directed = FALSE)
      }

      # Barabasi-Albert networks
      if (method_networks == "sf") {
        g = igraph::sample_pa(cluster_size, directed = FALSE)
      }

      # Converts the igraph object into an adjacency matrix
      adj <- as.matrix(igraph::as_adjacency_matrix(g))
      
      # Places this cluster's adjacency matrix into the appropriate block 
      # of the full matrix
      comba[idx, idx] <- adj
    }
    
    else {
      stop(paste0("Unknown method_networks: '", method_networks, "'. Must be one of 'ergm', 'er', 'sf'."))
    }
  }
  
  # Return result
  return(comba)
}
