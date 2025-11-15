#' @title
#' Identification of Partitions of the Network Causal Tree
#'
#' @description
#' Identifies the partitions of a Network Causal Tree
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  X N x M matrix, Observed Covariate Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize Minimum number of observations for each level of the joint intervention
#' to be required in the leafs
#' @param depth Depth of the tree
#'
#' @return A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the  Objective Function in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' valyes of the Xs to identify the given partition, the column TERMINAL reports the
#' 'role' of the node - parent or leaf -
#'
identify_partitions_nct <- function(method, alpha, beta, gamma,
                                    delta, depth, minsize,
                                    N, W, G, Y, X, p,
                                    Ne, Ne_list, population_effects) {
  
  # Initialize
  data_tree <- data.frame(idunit = 1:N, W = W, G = G, Y = Y)
  data_tree <- cbind(data_tree, as.data.frame(X))
  do_splits <- TRUE
  total_variance = NULL
  
  # Create output dataset
  tree_info <- data.frame(NODE = 1, OF = 0, NOBS = nrow(data_tree),
                          FILTER = NA, TERMINAL = "SPLIT",
                          stringsAsFactors = FALSE)
  
  # Recursively split the sample until one stopping criterion is met
  while (do_splits) {
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")
    # loop through every node
    for (j in to_calculate) {
      # extract that node's subset of data
      if (!is.na(tree_info[j, "FILTER"])) {
        texts = tree_info[j, "FILTER"]
        this_data <- subset(data_tree, eval(parse(text = texts)))
      } else {this_data <- data_tree}
      
      # number of nodes
      nleafs = nrow(tree_info)
      
      # find the best split for this node
      splitting <- compute_OF_Split( method = method, alpha = alpha, beta = beta,
                                     gamma = gamma, delta = delta,
                                     N = nrow(this_data), W = this_data$W,
                                     G = this_data$G, Y = this_data$Y,
                                     X = this_data[, grepl("^X\\.", names(this_data)), drop = FALSE],
                                     Ne = Ne[this_data$idunit], p = p[this_data$idunit],
                                     Ne_list = Ne_list[this_data$idunit],
                                     population_effects = population_effects,
                                     nleafs = nleafs, total_variance = total_variance)
      
      if (any(is.na(splitting))) {
        # node becomes a leaf
        split_here <- rep(FALSE, 2)
        
      } else {
        
        # extract split info
        maxof <- as.numeric(splitting[1])
        mn <- max(tree_info$NODE)
        
        # Paste splitting rules
        tmp_filter <- c(
          paste0(splitting[3], " >= ", as.numeric(splitting[2])),
          paste0(splitting[3], " < ", as.numeric(splitting[2]))
        )
      }
      
      
      # Check if the current splitting rule has already been used
      split_here  <- !sapply(tmp_filter,
                             FUN = function(x, y)
                               any(grepl(x, x = y)),
                             y = tree_info$FILTER)
      
      # Append splitting rules
      if (!is.na(tree_info[j, "FILTER"])) {
        tmp_filter  <- paste(tree_info[j, "FILTER"],
                             tmp_filter, sep = " & ") }
      
      # Check the number of observations in the current node
      tmp_nobs <- sapply(tmp_filter,
                         FUN = function(i, x){
                           nrow(subset(x = x, subset = eval(parse(text = i))))} ,
                         x = this_data)
      
      
      # STOP if insufficient minimum number of observations in the leafs
      if (any(as.numeric(table(this_data$W, this_data$G)) < minsize)) {
        split_here <- rep(FALSE, 2)
      }
      
      
      # STOP if exceeded depth
      depth_tree <- as.numeric(stringi::stri_count_regex(tree_info[j, "FILTER"], "X\\."))
      if (depth_tree >= depth & !is.na(depth_tree)) {
        split_here <- rep(FALSE, 2)
      }
      
      
      
      # Create children dataset
      children <- data.frame(NODE = c(mn + 1, mn + 2),
                             OF = c(rep(maxof,2)),
                             NOBS = tmp_nobs,
                             FILTER = tmp_filter,
                             TERMINAL = rep("SPLIT", 2),
                             row.names = NULL) [split_here,]
      
      if (method=="penalized") {
        
        variance_children = alpha*(Vartau1000(N = nrow(this_data), W = this_data$W,
                                              Y = this_data$Y, G = this_data$G,
                                              p = p[this_data$idunit], Ne = Ne[this_data$idunit],
                                              Ne_list = Ne_list[this_data$idunit])) +
          beta*(Vartau1101(N = nrow(this_data), W = this_data$W,
                           Y = this_data$Y, G = this_data$G,
                           p = p[this_data$idunit], Ne = Ne[this_data$idunit],
                           Ne_list = Ne_list[this_data$idunit])) +
          gamma*(Vartau1110(N = nrow(this_data), W = this_data$W,
                            Y = this_data$Y, G = this_data$G,
                            p = p[this_data$idunit], Ne = Ne[this_data$idunit],
                            Ne_list = Ne_list[this_data$idunit])) +
          delta*(Vartau0100(N = nrow(this_data), W = this_data$W,
                            Y = this_data$Y, G = this_data$G,
                            p = p[this_data$idunit], Ne = Ne[this_data$idunit],
                            Ne_list = Ne_list[this_data$idunit]))
        
        total_variance = c(total_variance, variance_children)
        
      }
      
      # update the tree
      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "SPLIT")
      tree_info <- rbind(tree_info, children)
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    }
    
  }
  
  return(tree = tree_info)
}
