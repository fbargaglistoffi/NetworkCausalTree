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
#' @param unit_ids Original network indices of the observations in W, G, and Y
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
compute_effects_nct <- function(output, nct_partition, N, W, G, Y, X,
                                Ne, Ne_list, p, minsize,
                                unit_ids = seq_len(N)) {
  
  if (is.null(X) || ncol(X) == 0) {
    X <- NULL
  } else {
    X <- as.data.frame(X)
    cn <- colnames(X)
    if (is.null(cn) || length(cn) != ncol(X) || !all(grepl("^X\\.", cn))) {
      colnames(X) <- paste0("X.", seq_len(ncol(X)))
    }
  }
  
  data_est <- data.frame(idunit = unit_ids, W = W, G = G, Y = Y)
  if (!is.null(X)) data_est <- cbind(data_est, X)
  if (!is.null(X)) {
    colnames(data_est) <- c("idunit", "W", "G", "Y", colnames(X))
  } else {
    colnames(data_est) <- c("idunit", "W", "G", "Y")
  }
  
  clean_filter_text = function(text) {
    text <- gsub("data_tree", "data_est", text)
    text <- gsub("\\bNA\\s*&\\s*", "", text)
    text <- gsub("\\bNA\\b", "", text)
    text <- gsub("&\\s*&", "&", text)
    text <- gsub("^&\\s*", "", text)
    text <- gsub("\\s*&$", "", text)
    trimws(text)
  }
  
  # If output equals to "estimation", then compute the estimated conditional average
  # treatment effects and their estimated variance, in all the partitions
  # identified by the tree
  
  if (output=="estimation") {

    NOBS_EST <- c(rep(0,nrow(nct_partition)))
    
    EFFTAU1000 = EFFTAU1101 = EFFTAU1110 = EFFTAU0100 = c(rep(0, nrow(nct_partition)))
    
    SETAU1000 = SETAU1101 = SETAU1110 = SETAU0100 = c(rep(0,nrow(nct_partition)))
    
    nct_partition <- cbind(nct_partition, NOBS_EST, EFFTAU1000, SETAU1000, EFFTAU1101,
                           SETAU1101, EFFTAU1110, SETAU1110, EFFTAU0100, SETAU0100)
    
    # Compute the effects
    
    nct_partition$NOBS_EST[1]<-N
    
    nct_partition$EFFTAU1000[1] <- EffTau1000(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU1101[1] <- EffTau1101(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU1110[1] <- EffTau1110(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU0100[1] <- EffTau0100(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    
    nct_partition$SETAU1000[1] <- sqrt(Vartau1000(N = nrow(data_est), W = data_est$W,
                                                  G = data_est$G, Y = data_est$Y, p = p,
                                                  Ne = Ne,
                                                  Ne_list = Ne_list,
                                                  unit_ids = data_est$idunit))
    nct_partition$SETAU1101[1] <- sqrt(Vartau1101(N = nrow(data_est), W = data_est$W,
                                                  G = data_est$G, Y = data_est$Y, p = p,
                                                  Ne = Ne,
                                                  Ne_list = Ne_list,
                                                  unit_ids = data_est$idunit))
    nct_partition$SETAU1110[1] <- sqrt(Vartau1110(N = nrow(data_est), W = data_est$W,
                                                  G = data_est$G, Y = data_est$Y, p = p,
                                                  Ne = Ne,
                                                  Ne_list = Ne_list,
                                                  unit_ids = data_est$idunit))
    nct_partition$SETAU0100[1] <- sqrt(Vartau0100(N = nrow(data_est), W = data_est$W,
                                                  G = data_est$G, Y = data_est$Y, p = p,
                                                  Ne = Ne,
                                                  Ne_list = Ne_list,
                                                  unit_ids = data_est$idunit))
    
    # loop over all the nodes except root
    if (nrow(nct_partition) > 1) {
      
      for (j in 2 : nrow(nct_partition)) {
        
        # Extract the filter condition
        if (!is.na(nct_partition[j, "FILTER"])) {
          
          texts <- clean_filter_text(nct_partition[j, "FILTER"])
        
          tmp <- with(data_est, eval(parse(text = texts)))

          if (sum(tmp) == 0) {
            nct_partition$NOBS_EST[j] <- 0
            nct_partition$EFFTAU1000[j] <- NA
            nct_partition$EFFTAU1101[j] <- NA
            nct_partition$EFFTAU1110[j] <- NA
            nct_partition$EFFTAU0100[j] <- NA
            nct_partition$SETAU1000[j] <- NA
            nct_partition$SETAU1101[j] <- NA
            nct_partition$SETAU1110[j] <- NA
            nct_partition$SETAU0100[j] <- NA
            next
          }
          
          this_data <- data_est[tmp, ]
          old_ids <- this_data$idunit
          
        } else {
          this_data <- data_est
          old_ids <- this_data$idunit
        }

        if (any(as.numeric(table(this_data$W, this_data$G)) < minsize)){
          warning('subpopulations not sufficiently represented')
        }

        nct_partition$NOBS_EST[j] <- nrow(this_data)
        
        # computed estimated effects
        nct_partition$EFFTAU1000[j] <- EffTau1000(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                  Y = this_data$Y, p = p,
                                                  Ne_list = Ne_list,
                                                  unit_ids = old_ids)
        nct_partition$EFFTAU1101[j] <- EffTau1101(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                  Y = this_data$Y, p = p,
                                                  Ne_list = Ne_list,
                                                  unit_ids = old_ids)
        nct_partition$EFFTAU1110[j] <- EffTau1110(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                  Y = this_data$Y, p = p,
                                                  Ne_list = Ne_list,
                                                  unit_ids = old_ids)
        nct_partition$EFFTAU0100[j] <- EffTau0100(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                  Y = this_data$Y, p = p,
                                                  Ne_list = Ne_list,
                                                  unit_ids = old_ids)
        
        
        # compute standard errors
        nct_partition$SETAU1000[j] <- sqrt(Vartau1000(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                      Y = this_data$Y, p = p,
                                                      Ne = Ne,
                                                      Ne_list = Ne_list,
                                                      unit_ids = old_ids))
        nct_partition$SETAU1101[j] <- sqrt(Vartau1101(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                      Y = this_data$Y, p = p,
                                                      Ne = Ne,
                                                      Ne_list = Ne_list,
                                                      unit_ids = old_ids))
        nct_partition$SETAU1110[j] <- sqrt(Vartau1110(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                      Y = this_data$Y, p = p,
                                                      Ne = Ne,
                                                      Ne_list = Ne_list,
                                                      unit_ids = old_ids))
        nct_partition$SETAU0100[j] <- sqrt(Vartau0100(N = nrow(this_data), W = this_data$W, G = this_data$G,
                                                      Y = this_data$Y, p = p,
                                                      Ne = Ne,
                                                      Ne_list = Ne_list,
                                                      unit_ids = old_ids))
      }
    }
    
    names(nct_partition)[names(nct_partition) == "SETAU1000"] <- "SE1000_EST"
    names(nct_partition)[names(nct_partition) == "SETAU1101"] <- "SE1101_EST"
    names(nct_partition)[names(nct_partition) == "SETAU1110"] <- "SE1110_EST"
    names(nct_partition)[names(nct_partition) == "SETAU0100"] <- "SE0100_EST"
    
  }
  
  # If output equals to "detection", then only compute the estimated conditional average
  # treatment effects  in all the partitions
  # identified by the tree
  
  if (output == "detection") {
    
    NOBS_EST = EFFTAU1000 = EFFTAU1101 = EFFTAU1110 = EFFTAU0100 = c(rep(0,nrow(nct_partition)))
    
    nct_partition <- cbind(nct_partition, NOBS_EST, EFFTAU1000, EFFTAU1101, EFFTAU1110, EFFTAU0100)
    
    nct_partition$NOBS_EST[1] <- N
    
    # Compute the effects
    
    nct_partition$EFFTAU1000[1] <- EffTau1000(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU1101[1] <- EffTau1101(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU1110[1] <- EffTau1110(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    nct_partition$EFFTAU0100[1] <- EffTau0100(N = nrow(data_est), W = data_est$W, G = data_est$G,
                                              Y = data_est$Y, p = p,
                                              Ne_list = Ne_list,
                                              unit_ids = data_est$idunit)
    
    if (nrow(nct_partition) > 1) {
      
      for (j in 2 : nrow(nct_partition)){
        
        
        if (!is.na(nct_partition[j, "FILTER"])) {
          
          texts <- clean_filter_text(nct_partition[j, "FILTER"])
          
          tmp <- with(data_est, eval(parse(text = texts)))

          if (sum(tmp) == 0) {
            nct_partition$NOBS_EST[j] <- 0
            nct_partition$EFFTAU1000[j] <- NA
            nct_partition$EFFTAU1101[j] <- NA
            nct_partition$EFFTAU1110[j] <- NA
            nct_partition$EFFTAU0100[j] <- NA
            next
          }
          
          this_data <- data_est[tmp, ]

          old_ids <- this_data$idunit
          
        } else {
          this_data <- data_est
          old_ids <- this_data$idunit
        }
        
        if (any(as.numeric(table(this_data$W, this_data$G)) < minsize)){
          warning('subpopulations not sufficiently represented')
        }
        
        nct_partition$NOBS_EST[j]<-nrow(this_data)
        nct_partition$EFFTAU1000[j] <- EffTau1000(N = nrow(this_data), W = this_data$W,
                                                  G = this_data$G, Y = this_data$Y,
                                                  p = p, Ne_list = Ne_list, unit_ids = old_ids)
        nct_partition$EFFTAU1101[j] <- EffTau1101(N = nrow(this_data), W = this_data$W,
                                                  G = this_data$G, Y = this_data$Y,
                                                  p = p, Ne_list = Ne_list, unit_ids = old_ids)
        nct_partition$EFFTAU1110[j] <- EffTau1110(N = nrow(this_data), W = this_data$W,
                                                  G = this_data$G, Y = this_data$Y,
                                                  p = p, Ne_list = Ne_list, unit_ids = old_ids)
        nct_partition$EFFTAU0100[j] <- EffTau0100(N = nrow(this_data), W = this_data$W,
                                                  G = this_data$G, Y = this_data$Y,
                                                  p = p, Ne_list = Ne_list, unit_ids = old_ids)
        
      } }
  }
  
  names(nct_partition)[names(nct_partition) == "NOBS"] <- "NOBS_TR"
  nct_partition$NOBS_TR <- as.numeric(nct_partition$NOBS_TR)
  nct_partition$TERMINAL <- as.character(nct_partition$TERMINAL)
  names(nct_partition)[names(nct_partition) == "EFFTAU1000"] <- "EFF1000_EST"
  names(nct_partition)[names(nct_partition) == "EFFTAU1101"] <- "EFF1101_EST"
  names(nct_partition)[names(nct_partition) == "EFFTAU1110"] <- "EFF1110_EST"
  names(nct_partition)[names(nct_partition) == "EFFTAU0100"] <- "EFF0100_EST"
  
  return(nct_partition)
}

