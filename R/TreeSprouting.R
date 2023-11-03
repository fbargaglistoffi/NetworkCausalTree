#' @title
#' Value of the Objective Function (OF)
#'
#' @description
#' Computes the measure of the Objective Function
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param total_variance - to be included if method = "penalized" - whole variance
#' @param nleafs - to be included if method = "penalized" - number of leafs
#'
#' @return A numeric value corresponding to the computed  Objective Function
#'
compute_OF_Value = function(method, alpha, beta, gamma, delta,
                   N, W, G, Y, p, Ne, Ne_list, population_effects,
                   total_variance, nleafs){

  # initialize
  inof=NULL

  # composite criterion
  if (method == "composite") {

    inof <- alpha * (((EffTau1000(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[1]) ^ 2) +
            beta * (((EffTau1101(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                      (population_effects[2]) ^ 2) +
            gamma * (((EffTau1110(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[3]) ^ 2) +
            delta * (((EffTau0100(N = N, W = W, G = G,Y = Y,p = p,Ne = Ne)) ^ 2) /
                       (population_effects[4]) ^ 2)
  }

  # penalized criterion
  if (method == "penalized") {

  inof <-  alpha * (((EffTau1000(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                    2 / nleafs * sum(c(total_variance, Vartau1000(N = N, W = W, Y = Y,G  = G,
                    p = p, Ne = Ne, Ne_list = Ne_list)))) +
           beta *  (((EffTau1101(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                    2 / nleafs *sum(c(total_variance, Vartau1101(N = N, W = W, Y = Y, G = G,
                    p = p, Ne = Ne, Ne_list = Ne_list)))) +
           gamma * (((EffTau1110( N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                   2 / nleafs * sum(c(total_variance,Vartau1110(N = N, W = W, Y = Y, G = G,
                   p = p, Ne = Ne, Ne_list = Ne_list))))+
           delta * (((EffTau0100(N = N, W = W, G = G, Y = Y, p = p, Ne = Ne)) ^ 2) -
                   2 / nleafs * sum(c(total_variance,Vartau0100(N = N, W = W, Y = Y, G = G,
                   p = p, Ne = Ne, Ne_list = Ne_list)))
    )
  }

  # singular criterion
  if (method == "singular") {

   inof <- alpha * ((EffTau1000( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2) +
           beta * ((EffTau1101( N = N, W = W, G = G,
                                Y = Y, p = p, Ne = Ne)) ^ 2) +
           gamma * ((EffTau1110( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2) +
           delta * ((EffTau0100( N = N, W = W, G = G,
                                 Y = Y, p = p, Ne = Ne)) ^ 2)
  }

  return(inof)
}



#-------------------------------------------------------------------------------

#' @title
#' Split Objective Function
#'
#' @description
#' Splits the sample where the  Objective Function is maximized
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param X N x K matrix, Observed Covariate Matrix
#' @param p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param total_variance - to be included if method = "penalized" - whole variance
#' @param nleafs - to be included if method = "penalized" - number of leafs
#'
#' @return A numeric vector made up by three elements: the first one identifies the
#' value of the  Objective Function is max , the second one reports the value of the variable that maximizes the  Objective Function,
#' the third one reports the corresponding variable
#'
compute_OF_Split = function(method, alpha, beta, gamma, delta,
                          N, W, G, Y, X, p, Ne, Ne_list,
                          population_effects, total_variance, nleafs){

  #initialize
  of<- c()
  name <-c()
  values<-c()

  # loop over all the predictors, identify all the unique splits for each variable,
  # and compute the corresponding value of OF

  for (j in 1:dim(X)[2]) {

    x = X[,j]
    splits <- sort(unique(x))
    valuesx<-c()
    ofx <- c()
    namesx <- c()

    # loop over all the possible splits

    for (i in seq_along(splits)) {

      sp <- splits[i]

      if (all(as.numeric(table(W[x >= sp], G[x >= sp]))>2) &
          all(as.numeric(table(W[x < sp], G[x < sp])) > 2)) {

        # Compute the corresponding OF
        ofx[i]<- 1/2*(compute_OF_Value(method = method, alpha = alpha, beta = beta,
                                       gamma = gamma,delta = delta, N = length(which(x < sp)),
                                       W = W[x < sp], G = G[x < sp],Y = Y[x < sp],
                                       Ne = Ne[x < sp], p = p[x < sp],
                                       population_effects = population_effects,
                                       Ne_list = Ne_list[x < sp],  nleafs = nleafs,
                                       total_variance = total_variance) +
                      compute_OF_Value(method = method, alpha = alpha, beta = beta,
                                       gamma = gamma,delta = delta, N = length(which(x >= sp)),
                                       W = W[x >= sp], G = G[x >= sp], Y = Y[x >= sp],
                                       Ne = Ne[x >= sp], p = p[x >= sp],
                                       population_effects = population_effects,
                                       Ne_list = Ne_list[x >= sp], nleafs = nleafs,
                                       total_variance = total_variance))
        }
      else {ofx[i] <- 0}
    }

    namex = rep(colnames(X)[j], length(unique(x)))
    valuesx = c(sort(unique(x)))

    # append all the computed values of the OF, all the variables that have defined the split
    # and their corresponding value.

    of = c(of, ofx)
    name = c(name, namex)
    values = c(values, valuesx)
  }

  # identify the variable and the exact splits which naximizes OF and the corresponding OF value

  if (all(is.na(of))) {
    ofres <- NA
    splitres <- NA
    varres <- NA
    } else {
    ofres <- max(na.omit(of))
    splitres <- values[which.max(of)]
    varres <- name[which.max(of)]
  }

    return(c(of = ofres , split = splitres , var=varres))

}


#-------------------------------------------------------------------------------

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
#' @param  X N x K matrix, Observed Covariate Matrix
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
  data_tree <- data.frame(idunit = c(1 : N) , W = W, G = G, Y = Y, X = X)
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

      splitting <- compute_OF_Split( method = method, alpha = alpha, beta = beta,
                                     gamma = gamma, delta = delta,
                                     N = nrow(this_data), W = this_data$W,
                                     G = this_data$G, Y = this_data$Y,
                                     X = this_data[, grepl("X.", names(this_data))],
                                     Ne = Ne[this_data$idunit], p = p[this_data$idunit],
                                     Ne_list = Ne_list[this_data$idunit],
                                     population_effects = population_effects,
                                     nleafs = nleafs, total_variance = total_variance)

      if (any(is.na(splitting))) {

        split_here <- rep(FALSE, 2)
        print('splits has stopped couse OF is all NA')

      } else {

        maxof <- as.numeric(splitting[1])
        mn <- max(tree_info$NODE)

      # Paste splitting rules
        tmp_filter <- c(paste("data_tree$",splitting[3], ">=","(" ,
                              as.numeric(splitting[2]),")",sep=""),
                        paste("data_tree$",splitting[3], "<", "(",
                              as.numeric(splitting[2]),")",sep=""))
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
        print('split has stopped for insufficient minsize')
      }


      # STOP if exceeded depth
      depth_tree <- as.numeric(stringi::stri_count_regex(tree_info[j, "FILTER"], "X."))
      if (depth_tree >= depth & !is.na(depth_tree)) {
        split_here <- rep(FALSE, 2)
        print('split has stopped for reached depth')
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

      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")
      tree_info <- rbind(tree_info, children)
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    }

  }

  return(tree = tree_info)
}


#-------------------------------------------------------------------------------

#' @title
#' Generation of the Network Causal Tree
#'
#' @description
#' Sprouts the network causal tree, eventually including a fraction of the initial -discovery-
#' sample.
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param  N Sample size
#' @param sampled_clusters Clusters assigned to the discovery set
#' @param m Total number of clusters
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  M N x 1 vector, Cluster Membership
#' @param  X N x K matrix, Observed Covariates Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  population_effects 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize minimum number of observaztions for each level of the joint intervention
#' to be required in the leafs
#' @param depth depth of the tree
#'
#' @return A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the OF in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' valyes of the Xs to identify the given partition,, the column TERMINAL reports the
#' 'role' of the node - parent or leaf -.
#'
sprout_nct = function(method, sampled_clusters,
                      m, alpha, beta, gamma, delta,
                      depth, minsize, N, W, G, Y, X, M, p, Ne,
                      population_effects, Ne_list){


  # Initialize
  data <- data.frame(idunit = c(1:N), W = W, G = G,
                     Y = Y, X = X, M = M)

  # Take only those observations that have been assigned to the discovery set
  datasample <- data[which(M %in% sampled_clusters), ]
  datasample <- datasample[order(datasample$idunit), ]
  sampleid <- unique(datasample$idunit)

  W = as.numeric(datasample$W)
  G = as.numeric(datasample$G)
  Y = as.numeric(datasample$Y)
  X = as.matrix(datasample[, -c(1 : 4, dim(datasample)[2])])
  colnames(X)=sub("X.","",colnames(X))

  # Identify partitions of the NCT on the discovery set
  tree_info <- identify_partitions_nct(method = method, alpha = alpha, beta = beta,
                                       gamma = gamma, delta = delta,
                                       N = length(sampleid), depth = depth,
                                       minsize = minsize,
                                       W = W, G = G,Y = Y, X = X,
                                       p = p[sampleid], Ne = Ne[sampleid],
                                       Ne_list = Ne_list[sampleid],
                                       population_effects = population_effects)

  return(tree_info)
}

#-------------------------------------------------------------------------------
#' @title
#' Computation of the Effects in all NCT partitions
#'
#' @description
#' Computes the estimates in all the partitions identified by the Network Causal Tree
#'
#' @param output Desired output of the analysis. if output = "detection" only point estimates
#' are computed, if output = "estimation" both estimated effects and variances are computed
#' @param nct_partition An NCT data frame
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  X N x K matrix, Observed Covariates Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param minsize minimum number of observations for each level of the joint intervention
#' to be required in the leafs
#'
#' @return A data.frame describing the obtained Network Causal Trees.
#' Each row represents a partition (of a specific tree) with 10/14 entries.
#' Columns summary:
#' - `OF`: value of the OF in the corresponding partition,
#' - `NOBS_TR`: number of training observations in the partition,
#' - `FILTER`: values of the covariates `X` that identify the partition,
#' - `NOBS_EST`: number of estimation observations in the partition,
#' - `EFF1000_EST`: estimated 1000 effects in the partitions,
#' - `EFF1101_EST`: estimated 1101 effects in the partitions,
#' - `EFF1110_EST`: estimated 1110 effects in the partitions,
#' - `EFF0100_EST`: estimated 0100 effects in the partitions.
#' Additional columns summary (only if output = "Estimation"):
#' - `SETAU1000_EST`: estimated std. error of the 1000 effect in the partition,
#' - `SETAU1101_EST`: estimated std. error of the 1101 effect in the partition,
#' - `SETAU1110_EST`: estimated std. error of the 1110 effect in the partition,
#' - `SETAU0100_EST`: estimated std. error of the 0100 effect in the partition.
#'
compute_effects_nct=function(output, nct_partition, N, W, G, Y, X,
                             Ne, Ne_list, p, minsize){

    # If output equals to "estimation", then compute the estimated conditional average
    # treatment effects and their estimated variance, in all the partitions
    # identified by the tree

    if (output=="estimation") {

      # Initialize

      data_est <- data.frame(idunit=c(1:N), W = W, G = G, Y = Y, X = X)

      NOBS_EST <- c(rep(0,nrow(nct_partition)))

      EFFTAU1000 = EFFTAU1101 = EFFTAU1110 = EFFTAU0100 = c(rep(0, nrow(nct_partition)))

      SETAU1000 = SETAU1101 = SETAU1110 = SETAU0100 = c(rep(0,nrow(nct_partition)))

      nct_partition <- cbind(nct_partition, NOBS_EST, EFFTAU1000, SETAU1000, EFFTAU1101,
                         SETAU1101, EFFTAU1110, SETAU1110, EFFTAU0100, SETAU0100)

      # Compute the effects

      nct_partition$NOBS_EST[1]<-N

      nct_partition$EFFTAU1000[1] <- EffTau1000(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU1101[1] <- EffTau1101(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU1110[1] <- EffTau1110(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU0100[1] <- EffTau0100(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])

      nct_partition$SETAU1000[1] <- sqrt(Vartau1000(N = nrow(data_est), W = data_est$W,
                                                G = data_est$G, Y = data_est$Y, p = p[data_est$idunit],
                                                Ne = Ne[data_est$idunit], Ne_list = Ne_list))
      nct_partition$SETAU1101[1] <- sqrt(Vartau1101(N = nrow(data_est), W = data_est$W,
                                                G = data_est$G, Y = data_est$Y, p = p[data_est$idunit],
                                                Ne = Ne[data_est$idunit], Ne_list = Ne_list))
      nct_partition$SETAU1110[1] <- sqrt(Vartau1110(N = nrow(data_est), W = data_est$W,
                                                G = data_est$G, Y = data_est$Y, p = p[data_est$idunit],
                                                Ne = Ne[data_est$idunit], Ne_list = Ne_list))
      nct_partition$SETAU0100[1] <- sqrt(Vartau0100(N = nrow(data_est), W = data_est$W,
                                                G = data_est$G, Y = data_est$Y, p = p[data_est$idunit],
                                                Ne = Ne[data_est$idunit], Ne_list = Ne_list))

      if (nrow(nct_partition) > 1) {

        for (j in 2 : nrow(nct_partition)) {

          texts = gsub(pattern = "data_tree", replacement = "data_est",nct_partition[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text = texts)))

          if (any(as.numeric(table(this_data$W, this_data$G)) < minsize / 2)) {
            print('subpopulations not sufficiently represented in some nodes of the Estimation Set')
          }


          Ne_listsub = Ne_list[this_data$idunit]
          nct_partition$NOBS_EST[j] <- nrow(this_data)
          nct_partition$EFFTAU1000[j] <- EffTau1000(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                Y = this_data$Y, p = p[this_data$idunit],
                                                Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU1101[j] <- EffTau1101(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                Y = this_data$Y, p = p[this_data$idunit],
                                                Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU1110[j] <- EffTau1110(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                Y = this_data$Y, p = p[this_data$idunit],
                                                Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU0100[j] <- EffTau0100(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                Y = this_data$Y, p = p[this_data$idunit],
                                                Ne = Ne[this_data$idunit])


         nct_partition$SETAU1000[j] <- sqrt(Vartau1000(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                   Y = this_data$Y, p = p[this_data$idunit],
                                                   Ne = Ne[this_data$idunit],
                                                   Ne_list = Ne_listsub))
         nct_partition$SETAU1101[j] <- sqrt(Vartau1101(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                   Y = this_data$Y, p = p[this_data$idunit],
                                                   Ne = Ne[this_data$idunit],
                                                   Ne_list = Ne_listsub))
         nct_partition$SETAU1110[j] <- sqrt(Vartau1110(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                   Y = this_data$Y, p = p[this_data$idunit],
                                                   Ne = Ne[this_data$idunit],
                                                   Ne_list = Ne_listsub))
         nct_partition$SETAU0100[j] <- sqrt(Vartau0100(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                   Y = this_data$Y, p = p[this_data$idunit],
                                                   Ne = Ne[this_data$idunit],
                                                   Ne_list = Ne_listsub))
        }
      }

      colnames(nct_partition)<-c("OF", "FILTER", "TERMINAL", "NOBS_TR", "NOBS_EST",
                             "EFF1000_EST", "SE1000_EST", "EFF1101_EST", "SE1101_EST",
                             "EFF1110_EST", "SE1110_EST", "EFF0100_EST", "SE0100_EST")

    }

  # If output equals to "detection", then only compute the estimated conditional average
  # treatment effects  in all the partitions
  # identified by the tree

    if (output == "detection") {

      # Initialize

      data_est <- data.frame(idunit = c(1:N), W = W, G = G, Y = Y, X = X)

      NOBS_EST = EFFTAU1000 = EFFTAU1101 = EFFTAU1110 = EFFTAU0100 = c(rep(0,nrow(nct_partition)))

      nct_partition <- cbind(nct_partition, NOBS_EST, EFFTAU1000, EFFTAU1101, EFFTAU1110, EFFTAU0100)

      nct_partition$NOBS_EST[1] <- N

      # Compute the effects

      nct_partition$EFFTAU1000[1] <- EffTau1000(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU1101[1] <- EffTau1101(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU1110[1] <- EffTau1110(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])
      nct_partition$EFFTAU0100[1] <- EffTau0100(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                            Y = data_est$Y, p = p[data_est$idunit],
                                            Ne = Ne[data_est$idunit])

      if (nrow(nct_partition) > 1) {

        for (j in 2 : nrow(nct_partition)){


          texts = gsub(pattern = "data_tree", replacement = "data_est",nct_partition[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text = texts)))

          if (any(as.numeric(table(this_data$W, this_data$G)) < 3)){
            warning('subpopulations not sufficiently represented')
          }

          nct_partition$NOBS_EST[j]<-nrow(this_data)
          nct_partition$EFFTAU1000[j] <- EffTau1000(N = nrow(this_data), W = this_data$W,
                                                G = this_data$G, Y = this_data$Y,
                                                p = p[this_data$idunit], Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU1101[j] <- EffTau1101(N = nrow(this_data), W = this_data$W,
                                                G = this_data$G, Y = this_data$Y,
                                                p = p[this_data$idunit], Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU1110[j] <- EffTau1110(N = nrow(this_data), W = this_data$W,
                                                G = this_data$G, Y = this_data$Y,
                                                p = p[this_data$idunit], Ne = Ne[this_data$idunit])
          nct_partition$EFFTAU0100[j] <- EffTau0100(N = nrow(this_data), W = this_data$W,
                                                G = this_data$G, Y = this_data$Y,
                                                p = p[this_data$idunit], Ne = Ne[this_data$idunit])

        } }

      colnames(nct_partition)<-c("OF", "FILTER", "TERMINAL", "NOBS_TR", "NOBS_EST",
                             "EFF1000_EST", "EFF1101_EST", "EFF1110_EST", "EFF0100_EST")

    }

  return(nct_partition)
}

