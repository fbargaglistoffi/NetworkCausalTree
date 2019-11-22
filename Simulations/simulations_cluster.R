###########################################################
####   SIMULATIONS FOR CAUSAL TREE WITH INTERFERENCE   ####
###########################################################

rm(list=ls())
#setwd("G:\\Il mio Drive\\Research\\Networks\\Draft Costanza\\networks\\Functions")
setwd("/home/falco.bargaglistoffi/Desktop/R_files/network_tree/")

# Upload Libraries for Parellel Computing
library(doParallel)
library(foreach)
library(randomizr)
library(MLmetrics)

###########################
## Initialize Parameters ##
###########################

# N: number of data points
# n_cov: number of variables 
# prob: treatment assignment probability
# rho: correlation within the covariates
# seq: Effect Size Magnitudes
# nsim: number of datasets created
# m: group indicator for clusters
N = 1500
n_cov = 10
prob = 0.4 
mu = rep(0, n_cov)
rho = 0
seq <- seq(0.1, 10.1, 1)
nsim = 100
m =  15

# Set up Parallel Computation
#setup parallel backend to use many processors
cl <- makeCluster(20)
registerDoParallel(cl)

# This function takes an arbitrary number of lists x all of which much have the same structure    
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

##################################
##    2 RULES &  1 MAIN EFFECT  ##
##################################

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    se.tau.1000 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0] <- i
      tau1000[x1==1 & x2==1] <- -i
      
      ## Generate Treatment Effects
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g)
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(1,0,0,0), method = "singular", # composite
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         se.est1000_1, se.est1000_2)
    }
    # Return the values
    list(correct.rules, tau.1000, se.tau.1000)
  }
})


correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
se_tau_1000 <- na.omit(matrix[[3]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))

nctree <- cbind(avg_correct_rules, mse_tau_1000, bias_tau_1000)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", 
                      "bias_tau_1000")
write.csv(nctree, file = "one_netcausaltree_main_effect.csv")


##################################
##    2 RULES & MAIN EFFECTS    ##
##################################

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    tau.1101 <- data.frame()
    se.tau.1000 <- data.frame()
    se.tau.1101 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0] <- i
      tau1000[x1==1 & x2==1] <- -i
      
      # Tau 1101
      tau1101 <- rep(0, n)
      tau1101[x1==0 & x2==0] <- i
      tau1101[x1==1 & x2==1] <- -i
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1101_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
          se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est1101_1 <- 0
          effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1101_1 <- 0
          se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est1101_1 <- 0
        effects.est1101_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est1101_1 <- 0
        se.est1101_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    # Return the values
    list(correct.rules, tau.1000, tau.1101, se.tau.1000, se.tau.1101)
  }
})


correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
tau_1101 <- na.omit(matrix[[3]])
se_tau_1000 <- na.omit(matrix[[4]])
se_tau_1101 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
tau_1101[tau_1101==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
mse_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))
bias_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V11) - 10.1) , na.rm = TRUE ))


nctree <- cbind(avg_correct_rules, mse_tau_1000, mse_tau_1101, bias_tau_1000, bias_tau_1101)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", "mse_tau_1101",
                      "bias_tau_1000", "bias_tau_1101")
write.csv(nctree, file = "netcausaltree_main_effects.csv")


##################################
##  2 RULES & SPILL, EFFECTS    ##
##################################

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    eta.1110 <- data.frame()
    eta.0100 <- data.frame()
    se.eta.1110 <- data.frame()
    se.eta.0100 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules eta. Direct effect=i over the population
      
      # eta1110
      eta1110 <- rep(0, n)
      eta1110[x1==0 & x2==0] <- i
      eta1110[x1==1 & x2==1] <- -i
      
      # eta 0100
      eta0100 <- rep(0, n)
      eta0100[x1==0 & x2==0] <- i
      eta0100[x1==1 & x2==1] <- -i
      
      ## Generate Treatment Effects
      y00 <- rnorm(n)
      y01 <- y00 + eta0100
      y10 <- rnorm(n)
      y11 <- y10 + eta1110
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0,0,0.5,0.5), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1110_1 <- SNCT$EFF1110_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1110_2 <- SNCT$EFF1110_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est0100_1 <- SNCT$EFF0100_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est0100_2 <- SNCT$EFF0100_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1110_1 <- SNCT$SE1110_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1110_2 <- SNCT$SE1110_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est0100_1 <- SNCT$SE0100_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est0100_2 <- SNCT$SE0100_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1110_1 <- SNCT$EFF1110_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1110_2 <- 0
          effects.est0100_1 <- SNCT$EFF0100_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est0100_2 <- 0
          se.est1110_1 <- SNCT$SE1110_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1110_2 <- 0
          se.est0100_1 <- SNCT$SE0100_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est0100_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1110_1 <- 0
          effects.est1110_2 <- SNCT$EFF1110_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est0100_1 <- 0
          effects.est0100_2 <- SNCT$EFF0100_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1110_1 <- 0
          se.est1110_2 <- SNCT$SE1110_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est0100_1 <- 0
          se.est0100_2 <- SNCT$SE0100_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1110_1 <- 0
        effects.est1110_2 <- 0
        effects.est0100_1 <- 0
        effects.est0100_2 <- 0
        se.est1110_1 <- 0
        se.est1110_2 <- 0
        se.est0100_1 <- 0
        se.est0100_2 <- 0
      }
      
      eta.1110[j*2-1, which(seq==i)] <- effects.est1110_1
      eta.1110[j*2, which(seq==i)] <- effects.est1110_2
      se.eta.1110[j*2-1, which(seq==i)] <-  se.est1110_1
      se.eta.1110[j*2, which(seq==i)] <-  se.est1110_2
      eta.0100[j*2-1, which(seq==i)] <- effects.est0100_1
      eta.0100[j*2, which(seq==i)] <- effects.est0100_2
      se.eta.0100[j*2-1, which(seq==i)] <-  se.est0100_1
      se.eta.0100[j*2, which(seq==i)] <-  se.est0100_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1110_1, effects.est1110_2,
         effects.est0100_1, effects.est0100_2,
         se.est1110_1, se.est1110_2,
         se.est0100_1, se.est0100_2)
    }
    # Return the values
    list(correct.rules, eta.1110, eta.0100, se.eta.1110, se.eta.0100)
  }
})


correct_rules <- na.omit(matrix[[1]])
eta_1110 <- na.omit(matrix[[2]])
eta_0100 <- na.omit(matrix[[3]])
se_eta_1110 <- na.omit(matrix[[4]])
se_eta_0100 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
eta_1110[eta_1110==0] <- NA
eta_0100[eta_0100==0] <- NA
mse_eta_1110 <- rbind(mean( (abs(eta_1110$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_1110$V11) - 10.1)^2 , na.rm = TRUE ))
mse_eta_0100 <- rbind(mean( (abs(eta_0100$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(eta_0100$V11) - 10.1)^2 , na.rm = TRUE ))
bias_eta_1110 <- rbind(mean( (abs(eta_1110$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(eta_1110$V11) - 10.1) , na.rm = TRUE ))
bias_eta_0100 <- rbind(mean( (abs(eta_0100$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(eta_0100$V11) - 10.1) , na.rm = TRUE ))


nctree_spillover <- cbind(avg_correct_rules, mse_eta_1110, mse_eta_0100, bias_eta_1110, bias_eta_0100)
colnames(nctree_spillover) <- c("correct_rules",
                                "mse_eta_1110", "mse_eta_0100",
                                "bias_eta_1110", "bias_eta_0100")
write.csv(nctree_spillover, file = "netcausaltree_spillover.csv")

##################################
##  3 RULES & SPILL, EFFECTS    ##
##################################

##################################
##    OVERLAP H.D. VARIABLES   ##
##################################

# In this design one rule has the same HDVs while the other had different ones.


set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    tau.1101 <- data.frame()
    se.tau.1000 <- data.frame()
    se.tau.1101 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0] <- i
      tau1000[x1==1 & x2==1] <- -i # different HDVs
      
      # Tau 1101
      tau1101 <- rep(0, n)
      tau1101[x1==0 & x2==0] <- i
      tau1101[x1==1 & x3==1] <- -i # different HDVs
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" | # change here
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                      rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          effects.est1101_2 <- 0 
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
          se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          se.est1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est1101_1 <- 0
          effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1101_1 <- 0
          se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est1101_1 <- 0
        effects.est1101_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est1101_1 <- 0
        se.est1101_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    # Return the values
    list(correct.rules, tau.1000, tau.1101, se.tau.1000, se.tau.1101)
  }
})
stopCluster()

correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
tau_1101 <- na.omit(matrix[[3]])
se_tau_1000 <- na.omit(matrix[[4]])
se_tau_1101 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
tau_1101[tau_1101==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
mse_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))
bias_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V11) - 10.1) , na.rm = TRUE ))


nctree <- cbind(avg_correct_rules, mse_tau_1000, mse_tau_1101, bias_tau_1000, bias_tau_1101)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", "mse_tau_1101",
                      "bias_tau_1000", "bias_tau_1101")
write.csv(nctree, file = "netcausaltree_main_effects_overlap.csv")

##################################
##    3 RULES & MAIN EFFECTS    ##
##################################

##################################
##     DIFFERENT EFFECT SIZES   ##
##################################

##################################
##    OVERLAP H.D. VARIABLES   ##
##################################

# In this design one rule has the same HDVs while the other had different ones.
# In addition one of the two effects is double of the other.

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    tau.1101 <- data.frame()
    se.tau.1000 <- data.frame()
    se.tau.1101 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0] <- i
      tau1000[x1==1 & x2==1] <- -i # different HDVs
      
      # Tau 1101
      tau1101 <- rep(0, n)
      tau1101[x1==0 & x2==0] <- i*2
      tau1101[x1==1 & x3==1] <- -i*2 # different HDVs
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" | # change here
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                      rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          effects.est1101_2 <- 0 
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
          se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          se.est1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est1101_1 <- 0
          effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1101_1 <- 0
          se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est1101_1 <- 0
        effects.est1101_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est1101_1 <- 0
        se.est1101_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    # Return the values
    list(correct.rules, tau.1000, tau.1101, se.tau.1000, se.tau.1101)
  }
})
stopCluster()

correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
tau_1101 <- na.omit(matrix[[3]])
se_tau_1000 <- na.omit(matrix[[4]])
se_tau_1101 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
tau_1101[tau_1101==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
mse_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))
bias_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V11) - 10.1) , na.rm = TRUE ))


nctree <- cbind(avg_correct_rules, mse_tau_1000, mse_tau_1101, bias_tau_1000, bias_tau_1101)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", "mse_tau_1101",
                      "bias_tau_1000", "bias_tau_1101")
write.csv(nctree, file = "netcausaltree_main_effects_overlap_diff_size.csv")

##################################
##    2 RULES & MAIN EFFECTS    ##
##################################

##################################
##      DEPTH OF THE RULES 3    ##
##################################

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    tau.1101 <- data.frame()
    se.tau.1000 <- data.frame()
    se.tau.1101 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0 & x3==0] <- i
      tau1000[x1==1 & x2==1 & x3==1] <- -i
      
      # Tau 1101
      tau1101 <- rep(0, n)
      tau1101[x1==0 & x2==0 & x3==0] <- i
      tau1101[x1==1 & x2==1 & x3==1] <- -i
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 3,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Introduce new combinations
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1 " |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1101_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
          se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est1101_1 <- 0
          effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1101_1 <- 0
          se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est1101_1 <- 0
        effects.est1101_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est1101_1 <- 0
        se.est1101_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    # Return the values
    list(correct.rules, tau.1000, tau.1101, se.tau.1000, se.tau.1101)
  }
})


correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
tau_1101 <- na.omit(matrix[[3]])
se_tau_1000 <- na.omit(matrix[[4]])
se_tau_1101 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
tau_1101[tau_1101==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
mse_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))
bias_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V11) - 10.1) , na.rm = TRUE ))


nctree <- cbind(avg_correct_rules, mse_tau_1000, mse_tau_1101, bias_tau_1000, bias_tau_1101)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", "mse_tau_1101",
                      "bias_tau_1000", "bias_tau_1101")
write.csv(nctree, file = "netcausaltree_main_effects_depth.csv")

##################################
##    2 RULES & MAIN EFFECTS    ##
##################################

##################################
##     CORRELATED REGRESSORS    ##
##################################

###########################
## Initialize Parameters ##
###########################


# rho: correlation within the covariates
rho = 0.7

set.seed(2019)
system.time({
  matrix <- foreach(j = 1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Load Packages and Functions
    source("speedfunctions.R")
    library(MASS)
    library(assertthat)
    library(dplyr)
    library(purrr)
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(stringi)
    library(gtools)
    library(tidyverse)
    
    ## Initialize Matrices
    correct.rules <- data.frame()
    tau.1000 <- data.frame()
    tau.1101 <- data.frame()
    se.tau.1000 <- data.frame()
    se.tau.1101 <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Cluster construction
      M <- randomizr::complete_ra(N = N, num_arms = m)
      levels(M) <- c(1:m)
      Mg <- as.numeric(table(M)) #groups size
      
      #Generate treatment
      p <- runif(m, min = prob, max = prob) #m dimensioned vector identifying the assignment prob. in each group
      
      # Assign individual assignment prob
      prt = c()
      for(k in 1:m){
        prt[which(M==k)] <- p[k]
      }
      
      # Randomly assign unit to treatment arms
      treat = c()
      for(k in 1:N){
        treat[k] <- rbinom(1, 1, prob=prt[k])
      }
      
      # Generate Variables
      Sigma = matrix(rho, nrow = n_cov, ncol = n_cov) + diag(n_cov)*(1-rho)
      rawvars = mvrnorm(n=N, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      
      # Generate outcome  variable
      outcome <- round(rnorm(N, mean = 20, sd = sqrt(10)), 2)
      
      # Adjacency matrix (generate Erdos-Reni within clusters)
      pl <- runif(m, min = 0.005, max = 0.01) 
      adiac_matrix <- matrix(0, N, N)
      for (k in 1:nrow(adiac_matrix)){
        for (q in 1:ncol(adiac_matrix)){
          if(k != q & M[k] == M[q]){
            adiac_matrix[k,q] <- rbinom(1, 1, prob = pl[M[k]])
          }  
        }
      }
      
      adiac_matrix[lower.tri(adiac_matrix)] <- t(adiac_matrix)[lower.tri(adiac_matrix)]
      neigh <- rowSums(adiac_matrix)
      
      net <- graph_from_adjacency_matrix(adiac_matrix, mode = "undirected") 
      
      #Compute number of treated neighbors and consequently G_i
      num_tr_neigh <- as.vector(adiac_matrix %*% treat) 
      neightreat <- rep(1, N) #G_i
      neightreat[num_tr_neigh==0] <- 0
      
      # Pass to the standard notation (exclude isoletad nodes)
      w <- treat[which(neigh != 0)]
      g <- neightreat[which(neigh != 0)]
      y <- outcome[which(neigh != 0)]
      M <- M[which(neigh != 0)]
      X <- X[which(neigh != 0),]
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      probT <- prt[which(neigh != 0)]
      NeighNum <- neigh[which(neigh != 0)] #degree
      n <- length(which(neigh != 0))
      adiac_matrix <- adiac_matrix[which(neigh != 0), which(neigh != 0)]
      
      
      ###################################################
      ## Generate Causal Rules Tau. Direct effect=i over the population
      
      # Tau1000
      tau1000 <- rep(0, n)
      tau1000[x1==0 & x2==0] <- i
      tau1000[x1==1 & x2==1] <- -i
      
      # Tau 1101
      tau1101 <- rep(0, n)
      tau1101[x1==0 & x2==0] <- i
      tau1101[x1==1 & x2==1] <- -i
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 7, mest = 8,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 30,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      if (correct==2) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                      rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                      rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1000_2 <- 0
          effects.est1101_1 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.est1101_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1000_2 <- 0
          se.est1101_1 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.est1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.est1101_1 <- 0
          effects.est1101_2 <- SNCT$EFF1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.est1101_1 <- 0
          se.est1101_2 <- SNCT$SE1101_EST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est1101_1 <- 0
        effects.est1101_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est1101_1 <- 0
        se.est1101_2 <- 0
      }
      
      tau.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    # Return the values
    list(correct.rules, tau.1000, tau.1101, se.tau.1000, se.tau.1101)
  }
})


correct_rules <- na.omit(matrix[[1]])
tau_1000 <- na.omit(matrix[[2]])
tau_1101 <- na.omit(matrix[[3]])
se_tau_1000 <- na.omit(matrix[[4]])
se_tau_1101 <- na.omit(matrix[[5]])

avg_correct_rules <- colMeans(correct_rules)
tau_1000[tau_1000==0] <- NA
tau_1101[tau_1101==0] <- NA
mse_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1000$V11) - 10.1)^2 , na.rm = TRUE ))
mse_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V2) - 1.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V3) - 2.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V4) - 3.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V5) - 4.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V6) - 5.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V7) - 6.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V8) - 7.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V9) - 8.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V10) - 9.1)^2 , na.rm = TRUE ),
                      mean( (abs(tau_1101$V11) - 10.1)^2 , na.rm = TRUE ))
bias_tau_1000 <- rbind(mean( (abs(tau_1000$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1000$V11) - 10.1) , na.rm = TRUE ))
bias_tau_1101 <- rbind(mean( (abs(tau_1101$V1) - 0.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V2) - 1.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V3) - 2.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V4) - 3.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V5) - 4.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V6) - 5.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V7) - 6.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V8) - 7.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V9) - 8.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V10) - 9.1) , na.rm = TRUE ),
                       mean( (abs(tau_1101$V11) - 10.1) , na.rm = TRUE ))


nctree <- cbind(avg_correct_rules, mse_tau_1000, mse_tau_1101, bias_tau_1000, bias_tau_1101)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_1000", "mse_tau_1101",
                      "bias_tau_1000", "bias_tau_1101")
write.csv(nctree, file = "netcausaltree_main_effects_correlated.csv")

stopCluster()

