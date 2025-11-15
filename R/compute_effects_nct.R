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
#' @param  X N x M matrix, Observed Covariates Matrix
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
    
    # loop over all the nodes except root
    if (nrow(nct_partition) > 1) {
      
      for (j in 2 : nrow(nct_partition)) {
        
        # extract the filter condition
        if (!is.na(nct_partition[j, "FILTER"])) {
          texts = gsub("data_tree", "data_est", nct_partition[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text = texts)))
        } else {
          this_data <- data_est
        }
        
        # check subgroup size
        if (any(as.numeric(table(this_data$W, this_data$G)) < minsize / 2)) {
        }
        
        
        # extract neighbor info
        Ne_listsub = Ne_list[this_data$idunit]
        nct_partition$NOBS_EST[j] <- nrow(this_data)
        
        # computed estimated effects
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
        
        
        # compute standard errors
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
    
    # rename column to NOBS_TR
    names(nct_partition)[names(nct_partition) == "NOBS"] <- "NOBS_TR"
    
    # expected final column names
    expected_names <- c(
      "NODE","OF","NOBS_TR","FILTER","TERMINAL",
      "NOBS_EST",
      "EFF1000_EST","SE1000_EST","EFF1101_EST","SE1101_EST",
      "EFF1110_EST","SE1110_EST","EFF0100_EST","SE0100_EST"
    )
    
    # only assign up to the number of existing columns
    colnames(nct_partition)[seq_len(min(ncol(nct_partition), length(expected_names)))] <- 
      expected_names[seq_len(min(ncol(nct_partition), length(expected_names)))]
    
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
        
        
        if (!is.na(nct_partition[j, "FILTER"])) {
          texts = gsub("data_tree", "data_est", nct_partition[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text = texts)))
        } else {
          this_data <- data_est
        }
        
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
    names(nct_partition)[names(nct_partition) == "NOBS"] <- "NOBS_TR"
    
    expected_names <- c(
      "NODE","OF","NOBS_TR","FILTER","TERMINAL",
      "NOBS_EST","EFF1000_EST","EFF1101_EST","EFF1110_EST","EFF0100_EST"
    )
    
    colnames(nct_partition)[seq_len(min(ncol(nct_partition), length(expected_names)))] <- 
      expected_names[seq_len(min(ncol(nct_partition), length(expected_names)))]
    
  }
  
  return(nct_partition)
}

