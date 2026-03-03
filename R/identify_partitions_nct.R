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
  
  X <- as.data.frame(X)
  
  if (nrow(X) == 0) {
    X <- matrix(nrow = N, ncol = 0)
  }
  
  data_tree <- cbind(data_tree, X)
  
  do_splits <- TRUE
  total_variance = NULL
  
  # Create output dataset
  tree_info <- data.frame(NODE = 1, OF = 0, NOBS = nrow(data_tree),
                          FILTER = NA, TERMINAL = "SPLIT",
                          stringsAsFactors = FALSE)
  
  # Recursively split the sample until one stopping criterion is met
  while (do_splits) {
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")

    for (j in to_calculate) {

      if (!is.na(tree_info[j, "FILTER"])) {
        texts = tree_info[j, "FILTER"]
        this_data <- subset(data_tree, eval(parse(text = texts)))
      } else {this_data <- data_tree}

      nleafs = nrow(tree_info)

      Xnode <- this_data[, grepl("^X\\.", names(this_data)), drop = FALSE]

      if (ncol(Xnode) == 0) {
        tree_info[j, "TERMINAL"] <- "LEAF"
        next
      }

      splitting <- compute_OF_Split(
        method = method, alpha = alpha, beta = beta,
        gamma = gamma, delta = delta,
        W = this_data$W,
        G = this_data$G, Y = this_data$Y,
        X = Xnode,
        Ne = Ne[this_data$idunit], p = p[this_data$idunit],
        Ne_list = Ne_list[this_data$idunit],
        population_effects = population_effects,
        nleafs = nleafs, total_variance = total_variance
      )

      tmp_filter <- NULL
      split_here <- rep(FALSE, 2)
      
      if (any(is.na(splitting))) {
        tree_info[j, "TERMINAL"] <- "LEAF"
        
      } else {
        maxof <- as.numeric(splitting[1])
        mn <- max(tree_info$NODE)
        
        tmp_filter <- c(
          paste0(splitting[3], " >= ", as.numeric(splitting[2])),
          paste0(splitting[3], " < ", as.numeric(splitting[2]))
        )
        
        if (!is.na(tree_info[j, "FILTER"])) {
          tmp_filter  <- paste(tree_info[j, "FILTER"],
                               tmp_filter, sep = " & ") 
        }
        
        split_here  <- !sapply(tmp_filter,
                               FUN = function(x, y) x %in% y,
                               y = tree_info$FILTER)
        
        tmp_nobs <- sapply(tmp_filter,
                           FUN = function(i, x){
                             nrow(subset(x = x, subset = eval(parse(text = i))))} ,
                           x = this_data)
        
        for (child_idx in 1:2) {
          if (split_here[child_idx]) {
            child_data <- subset(this_data, eval(parse(text = tmp_filter[child_idx])))
            child_table <- table(child_data$W, child_data$G)
            
            if (any(as.numeric(child_table) < minsize)) {
              split_here[child_idx] <- FALSE
            }
          }
        }
        
        depth_tree <- as.numeric(stringi::stri_count_regex(tree_info[j, "FILTER"], "X\\."))
        if (depth_tree >= depth & !is.na(depth_tree)) {
          split_here <- rep(FALSE, 2)
        }
        
        if (any(split_here)) {
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
          
          tree_info <- rbind(tree_info, children)
          
          tree_info[j, "TERMINAL"] <- "PARENT"
        } else {
          tree_info[j, "TERMINAL"] <- "LEAF"
        }
      }
    }
    
    do_splits <- any(tree_info$TERMINAL == "SPLIT")
  }
  
  tree_info$TERMINAL[tree_info$TERMINAL == "PARENT"] <- "SPLIT"
  
  return(tree = tree_info)
}
