###########################################################
####   SIMULATIONS FOR CAUSAL TREE WITH INTERFERENCE   ####
###########################################################

rm(list=ls())
#setwd("D:\\Research\\Network Causal Tree\\Network-Causal-Tree\\Functions")
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
cl <- makeCluster(15)
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
    tau.est.1000 <- data.frame()
    se.tau.est.1000 <- data.frame()
    tau.test.1000 <- data.frame()
    se.tau.test.1000 <- data.frame()
    vi.x1 <- vi.x2 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
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
      SNCT <- SimNetworkCausalTrees(effweights = c(1,0,0,0), method = "singular", # composite
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 5, mest = 5, 
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 10,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3") | 
                       str_detect(SNCT$FILTER, "X.4") |
                       str_detect(SNCT$FILTER, "X.5") |
                       str_detect(SNCT$FILTER, "X.6") |
                       str_detect(SNCT$FILTER, "X.7") |
                       str_detect(SNCT$FILTER, "X.8") |
                       str_detect(SNCT$FILTER, "X.9") |
                       str_detect(SNCT$FILTER, "X.10")])
      
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      
      
      ## Extract values for the Test Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1000_2 <- 0
          se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1000_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1000_1 <- 0
          effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1000_1 <- 0
          se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1000_1 <- 0
        effects.test1000_2 <- 0
        se.test1000_1 <- 0
        se.test1000_2 <- 0
      }
      
      tau.test.1000[j*2-1, which(seq==i)] <- effects.test1000_1
      tau.test.1000[j*2, which(seq==i)] <- effects.test1000_2
      se.tau.test.1000[j*2-1, which(seq==i)] <-  se.test1000_1
      se.tau.test.1000[j*2, which(seq==i)] <-  se.test1000_2
      
      ## Remove Vectors
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         se.est1000_1, se.est1000_2)
    }
    
    ## Return the values
    
    list(correct.rules, tau.est.1000, se.tau.est.1000, tau.test.1000, se.tau.test.1000, vi.x1, vi.x2, vi.X)
  }
})

## Extract the Returned Values

correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
se_tau_est_1000 <- na.omit(matrix[[3]])
tau_test_1000 <- na.omit(matrix[[4]])
se_tau_test_1000 <- na.omit(matrix[[5]])
vi_x1 <- na.omit(matrix[[6]])
vi_x2 <- na.omit(matrix[[7]])
vi_X <- na.omit(matrix[[8]])


## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_X <- colMeans(vi_X)

# Exclude Rules that were not discovered

tau_est_1000[tau_est_1000==0] <- NA 
tau_test_1000[tau_test_1000==0] <- NA 

## Mean Squared Error (Estimation/Test Set)
mse_tau_est_1000 <- c()
mse_tau_test_1000 <- c()
for(i in seq){
  mse_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_tau_est_1000 <- c()
bias_tau_test_1000 <- c()
for(i in seq){
  bias_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1000 <- data.frame()
coverage_test_1000 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_est_1000)) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000[j,which(seq==i)] <- between(i, abs(tau_est_1000[j,which(seq==i)]) - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                  abs(tau_est_1000[j,which(seq==i)]) + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1000)) {
    if (is.na(tau_test_1000[j,which(seq==i)])==FALSE) {
      coverage_test_1000[j,which(seq==i)] <- between(i, abs(tau_test_1000[j,which(seq==i)]) - 1.96*se_tau_test_1000[j,which(seq==i)],
                                                    abs(tau_test_1000[j,which(seq==i)]) + 1.96*se_tau_test_1000[j,which(seq==i)])
    }
    else {
      coverage_test_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1000 <- colMeans(coverage_test_1000, na.rm = TRUE) 

## Create a Matrix for the Results

results_nctree <- cbind(avg_correct_rules,
                mse_tau_est_1000, bias_tau_est_1000,
                mse_tau_test_1000, bias_tau_test_1000,
                coverage_est_1000, coverage_test_1000, 
                vi_x1, vi_x2, vi_X)
colnames(results_nctree) <- c("correct_rules",
                      "mse_tau_est_1000", 
                      "bias_tau_est_1000",
                      "mse_tau_test_1000", 
                      "bias_tau_test_1000",
                      "coverage_est_1000",
                      "coverage_test_1000",
                      "vi_x1", "vi_x2", "vi_X")

## Save the Results

write.csv(results_nctree, file = "one_main_effect.csv")


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
    tau.est.1000 <- data.frame()
    tau.est.1101 <- data.frame()
    se.tau.est.1000 <- data.frame()
    se.tau.est.1101 <- data.frame()
    tau.test.1000 <- data.frame()
    tau.test.1101 <- data.frame()
    se.tau.test.1000 <- data.frame()
    se.tau.test.1101 <- data.frame()
    vi.x1 <- vi.x2 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
      # Compute number of treated neighbors and consequently G_i
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
      SNCT <- SimNetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 5, mest = 5,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 20,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3") | 
                                               str_detect(SNCT$FILTER, "X.4") |
                                               str_detect(SNCT$FILTER, "X.5") |
                                               str_detect(SNCT$FILTER, "X.6") |
                                               str_detect(SNCT$FILTER, "X.7") |
                                               str_detect(SNCT$FILTER, "X.8") |
                                               str_detect(SNCT$FILTER, "X.9") |
                                               str_detect(SNCT$FILTER, "X.10")])
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.est.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.est.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.est.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.est.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      ## Extract Values for the testimation Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1000_2 <- 0
          effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1101_2 <- 0
          se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1000_2 <- 0
          se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1000_1 <- 0
          effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.test1101_1 <- 0
          effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1000_1 <- 0
          se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1101_1 <- 0
          se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1000_1 <- 0
        effects.test1000_2 <- 0
        effects.test1101_1 <- 0
        effects.test1101_2 <- 0
        se.test1000_1 <- 0
        se.test1000_2 <- 0
        se.test1101_1 <- 0
        se.test1101_2 <- 0
      }
      
      tau.test.1000[j*2-1, which(seq==i)] <- effects.test1000_1
      tau.test.1000[j*2, which(seq==i)] <- effects.test1000_2
      se.tau.test.1000[j*2-1, which(seq==i)] <-  se.test1000_1
      se.tau.test.1000[j*2, which(seq==i)] <-  se.test1000_2
      tau.test.1101[j*2-1, which(seq==i)] <- effects.test1101_1
      tau.test.1101[j*2, which(seq==i)] <- effects.test1101_2
      se.tau.test.1101[j*2-1, which(seq==i)] <-  se.test1101_1
      se.tau.test.1101[j*2, which(seq==i)] <-  se.test1101_2
      
      ## Clear Memory
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    
    ## Return the values
    
    list(correct.rules, tau.est.1000, tau.est.1101, se.tau.est.1000, se.tau.est.1101,
         tau.test.1000, tau.test.1101, se.tau.test.1000, se.tau.test.1101, vi.x1, vi.x2, vi.X)
  }
})

## Extract the Results 

correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
tau_est_1101 <- na.omit(matrix[[3]])
se_tau_est_1000 <- na.omit(matrix[[4]])
se_tau_est_1101 <- na.omit(matrix[[5]])
tau_test_1000 <- na.omit(matrix[[6]])
tau_test_1101 <- na.omit(matrix[[7]])
se_tau_test_1000 <- na.omit(matrix[[8]])
se_tau_test_1101 <- na.omit(matrix[[9]])
vi_x1 <- na.omit(matrix[[10]])
vi_x2 <- na.omit(matrix[[11]])
vi_X <- na.omit(matrix[[12]])


## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered

tau_est_1000[tau_est_1000==0] <- NA
tau_est_1101[tau_est_1101==0] <- NA
tau_test_1000[tau_est_1000==0] <- NA
tau_test_1101[tau_est_1101==0] <- NA


## Mean Squared Error (Estimation/Test Set)

mse_tau_est_1000 <- c()
mse_tau_test_1000 <- c()
mse_tau_est_1101 <- c()
mse_tau_test_1101 <- c()

for(i in seq){
  mse_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_tau_est_1000 <- c()
bias_tau_test_1000 <- c()
bias_tau_est_1101 <- c()
bias_tau_test_1101 <- c()

for(i in seq){
  bias_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1000 <- data.frame()
coverage_test_1000 <- data.frame()
coverage_est_1101 <- data.frame()
coverage_test_1101 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_est_1000)) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000[j,which(seq==i)] <- between(i, abs(tau_est_1000[j,which(seq==i)]) - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                    abs(tau_est_1000[j,which(seq==i)]) + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1000)) {
    if (is.na(tau_test_1000[j,which(seq==i)])==FALSE) {
      coverage_test_1000[j,which(seq==i)] <- between(i, abs(tau_test_1000[j,which(seq==i)]) - 1.96*se_tau_test_1000[j,which(seq==i)],
                                                     abs(tau_test_1000[j,which(seq==i)]) + 1.96*se_tau_test_1000[j,which(seq==i)])
    }
    else {
      coverage_test_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1000 <- colMeans(coverage_test_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_est_1101)) {
    if (is.na(tau_est_1101[j,which(seq==i)])==FALSE) {
      coverage_est_1101[j,which(seq==i)] <- between(i, abs(tau_est_1101[j,which(seq==i)]) - 1.96*se_tau_est_1101[j,which(seq==i)],
                                                    abs(tau_est_1101[j,which(seq==i)]) + 1.96*se_tau_est_1101[j,which(seq==i)])
    }
    else {
      coverage_est_1101[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1101 <- colMeans(coverage_est_1101, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1101)) {
    if (is.na(tau_test_1101[j,which(seq==i)])==FALSE) {
      coverage_test_1101[j,which(seq==i)] <- between(i, abs(tau_test_1101[j,which(seq==i)]) - 1.96*se_tau_test_1101[j,which(seq==i)],
                                                     abs(tau_test_1101[j,which(seq==i)]) + 1.96*se_tau_test_1101[j,which(seq==i)])
    }
    else {
      coverage_test_1101[j,which(seq==i)] <- NA
    }
  }  
}

## Create a Matrix for the Results

nctree <- cbind(avg_correct_rules, mse_tau_est_1000, mse_tau_est_1101, bias_tau_est_1000, bias_tau_est_1101, coverage_est_1000, coverage_est_1101,
                mse_tau_test_1000, mse_tau_test_1101, bias_tau_test_1000, bias_tau_test_1101, coverage_test_1000, coverage_test_1101, vi_x1, vi_x2, vi_X)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_est_1000", "mse_tau_est_1101",
                      "bias_tau_est_1000", "bias_tau_est_1101",
                      "coverage_est_1000", "coverage_est_1101",
                      "mse_tau_test_1000", "mse_tau_test_1101",
                      "bias_tau_test_1000", "bias_tau_test_1101",
                      "coverage_test_1000", "coverage_test_1101",
                      "vi_x1", "vi_x2", "vi_X")
write.csv(nctree, file = "two_main_effects.csv")


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
    vi.x1 <- vi.x2 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
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
      SNCT <- SimNetworkCausalTrees(effweights = c(0,0,0.5,0.5), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = 5, mest = 5,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize = 20,
                                 n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3") | 
                                               str_detect(SNCT$FILTER, "X.4") |
                                               str_detect(SNCT$FILTER, "X.5") |
                                               str_detect(SNCT$FILTER, "X.6") |
                                               str_detect(SNCT$FILTER, "X.7") |
                                               str_detect(SNCT$FILTER, "X.8") |
                                               str_detect(SNCT$FILTER, "X.9") |
                                               str_detect(SNCT$FILTER, "X.10")])
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      eta.est.1110[j*2-1, which(seq==i)] <- effects.est1110_1
      eta.est.1110[j*2, which(seq==i)] <- effects.est1110_2
      se.eta.est.1110[j*2-1, which(seq==i)] <-  se.est1110_1
      se.eta.est.1110[j*2, which(seq==i)] <-  se.est1110_2
      eta.est.0100[j*2-1, which(seq==i)] <- effects.est0100_1
      eta.est.0100[j*2, which(seq==i)] <- effects.est0100_2
      se.eta.est.0100[j*2-1, which(seq==i)] <-  se.est0100_1
      se.eta.est.0100[j*2, which(seq==i)] <-  se.est0100_2
      
      ## Extract Values for the testimation Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1110_1 <- SNCT$EFF1110_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1110_2 <- SNCT$EFF1110_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.test0100_1 <- SNCT$EFF0100_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test0100_2 <- SNCT$EFF0100_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1110_1 <- SNCT$SE1110_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1110_2 <- SNCT$SE1110_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test0100_1 <- SNCT$SE0100_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test0100_2 <- SNCT$SE0100_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1110_1 <- SNCT$EFF1110_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1110_2 <- 0
          effects.test0100_1 <- SNCT$EFF0100_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test0100_2 <- 0
          se.test1110_1 <- SNCT$SE1110_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1110_2 <- 0
          se.test0100_1 <- SNCT$SE0100_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test0100_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1110_1 <- 0
          effects.test1110_2 <- SNCT$EFF1110_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.test0100_1 <- 0
          effects.test0100_2 <- SNCT$EFF0100_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1110_1 <- 0
          se.test1110_2 <- SNCT$SE1110_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test0100_1 <- 0
          se.test0100_2 <- SNCT$SE0100_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1110_1 <- 0
        effects.test1110_2 <- 0
        effects.test0100_1 <- 0
        effects.test0100_2 <- 0
        se.test1110_1 <- 0
        se.test1110_2 <- 0
        se.test0100_1 <- 0
        se.test0100_2 <- 0
      }
      
      eta.test.1110[j*2-1, which(seq==i)] <- effects.test1110_1
      eta.test.1110[j*2, which(seq==i)] <- effects.test1110_2
      se.eta.test.1110[j*2-1, which(seq==i)] <-  se.test1110_1
      se.eta.test.1110[j*2, which(seq==i)] <-  se.test1110_2
      eta.test.0100[j*2-1, which(seq==i)] <- effects.test0100_1
      eta.test.0100[j*2, which(seq==i)] <- effects.test0100_2
      se.eta.test.0100[j*2-1, which(seq==i)] <-  se.test0100_1
      se.eta.test.0100[j*2, which(seq==i)] <-  se.test0100_2
      
      ## Clear Memory
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1110_1, effects.est1110_2,
         effects.est0100_1, effects.est0100_2,
         se.est1110_1, se.est1110_2,
         se.est0100_1, se.est0100_2)
    }
    
    ## Return the values
    
    list(correct.rules, eta.est.1110, eta.est.0100, se.eta.est.1110, se.eta.est.0100,
         eta.test.1110, eta.test.0100, se.eta.test.1110, se.eta.test.0100, vi.x1, vi.x2, vi.X)
  }
})

## Extract the Results 

correct_rules <- na.omit(matrix[[1]])
eta_est_1110 <- na.omit(matrix[[2]])
eta_est_0100 <- na.omit(matrix[[3]])
se_eta_est_1110 <- na.omit(matrix[[4]])
se_eta_est_0100 <- na.omit(matrix[[5]])
eta_test_1110 <- na.omit(matrix[[6]])
eta_test_0100 <- na.omit(matrix[[7]])
se_eta_test_1110 <- na.omit(matrix[[8]])
se_eta_test_0100 <- na.omit(matrix[[9]])
vi_x1 <- na.omit(matrix[[10]])
vi_x2 <- na.omit(matrix[[11]])
vi_X <- na.omit(matrix[[12]])


## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered

eta_est_1110[eta_est_1110==0] <- NA
eta_est_0100[eta_est_0100==0] <- NA
eta_test_1110[eta_est_1110==0] <- NA
eta_test_0100[eta_est_0100==0] <- NA


## Mean Squared Error (Estimation/Test Set)

mse_eta_est_1110 <- c()
mse_eta_test_1110 <- c()
mse_eta_est_0100 <- c()
mse_eta_test_0100 <- c()

for(i in seq){
  mse_eta_est_1110[which(seq==i)] <- mean( (abs(eta_est_1110[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_eta_test_1110[which(seq==i)] <- mean( (abs(eta_test_1110[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_eta_est_0100[which(seq==i)] <- mean( (abs(eta_est_0100[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_eta_test_0100[which(seq==i)] <- mean( (abs(eta_test_0100[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_eta_est_1110 <- c()
bias_eta_test_1110 <- c()
bias_eta_est_0100 <- c()
bias_eta_test_0100 <- c()

for(i in seq){
  bias_eta_est_1110[which(seq==i)] <- mean( (abs(eta_est_1110[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_eta_test_1110[which(seq==i)] <- mean( (abs(eta_test_1110[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_eta_est_0100[which(seq==i)] <- mean( (abs(eta_est_0100[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_eta_test_0100[which(seq==i)] <- mean( (abs(eta_test_0100[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1110 <- data.frame()
coverage_test_1110 <- data.frame()
coverage_est_0100 <- data.frame()
coverage_test_0100 <- data.frame()

for(i in seq){
  for(j in 1:nrow(eta_est_1110)) {
    if (is.na(eta_est_1110[j,which(seq==i)])==FALSE) {
      coverage_est_1110[j,which(seq==i)] <- between(i, abs(eta_est_1110[j,which(seq==i)]) - 1.96*se_eta_est_1110[j,which(seq==i)],
                                                    abs(eta_est_1110[j,which(seq==i)]) + 1.96*se_eta_est_1110[j,which(seq==i)])
    }
    else {
      coverage_est_1110[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1110 <- colMeans(coverage_est_1110, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(eta_test_1110)) {
    if (is.na(eta_test_1110[j,which(seq==i)])==FALSE) {
      coverage_test_1110[j,which(seq==i)] <- between(i, abs(eta_test_1110[j,which(seq==i)]) - 1.96*se_eta_test_1110[j,which(seq==i)],
                                                     abs(eta_test_1110[j,which(seq==i)]) + 1.96*se_eta_test_1110[j,which(seq==i)])
    }
    else {
      coverage_test_1110[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1110 <- colMeans(coverage_test_1110, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(eta_est_0100)) {
    if (is.na(eta_est_0100[j,which(seq==i)])==FALSE) {
      coverage_est_0100[j,which(seq==i)] <- between(i, abs(eta_est_0100[j,which(seq==i)]) - 1.96*se_eta_est_0100[j,which(seq==i)],
                                                    abs(eta_est_0100[j,which(seq==i)]) + 1.96*se_eta_est_0100[j,which(seq==i)])
    }
    else {
      coverage_est_0100[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_0100 <- colMeans(coverage_est_0100, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(eta_test_0100)) {
    if (is.na(eta_test_0100[j,which(seq==i)])==FALSE) {
      coverage_test_0100[j,which(seq==i)] <- between(i, abs(eta_test_0100[j,which(seq==i)]) - 1.96*se_eta_test_0100[j,which(seq==i)],
                                                     abs(eta_test_0100[j,which(seq==i)]) + 1.96*se_eta_test_0100[j,which(seq==i)])
    }
    else {
      coverage_test_0100[j,which(seq==i)] <- NA
    }
  }  
}

## Create a Matrix for the Results

nctree <- cbind(avg_correct_rules, mse_eta_est_1110, mse_eta_est_0100, bias_eta_est_1110, bias_eta_est_0100, coverage_est_1110, coverage_est_0100,
                mse_eta_test_1110, mse_eta_test_0100, bias_eta_test_1110, bias_eta_test_0100, coverage_test_1110, coverage_test_0100, vi_x1, vi_x2, vi_X)
colnames(nctree) <- c("correct_rules",
                      "mse_eta_est_1110", "mse_eta_est_0100",
                      "bias_eta_est_1110", "bias_eta_est_0100",
                      "coverage_est_1110", "coverage_est_0100",
                      "mse_eta_test_1110", "mse_eta_test_0100",
                      "bias_eta_test_1110", "bias_eta_test_0100",
                      "coverage_test_1110", "coverage_test_0100",
                      "vi_x1", "vi_x2", "vi_X")
write.csv(nctree, file = "two_spillover_effects.csv")


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
    tau.est.1000 <- data.frame()
    tau.est.1101 <- data.frame()
    se.tau.est.1000 <- data.frame()
    se.tau.est.1101 <- data.frame()
    tau.test.1000 <- data.frame()
    tau.test.1101 <- data.frame()
    se.tau.test.1000 <- data.frame()
    se.tau.test.1101 <- data.frame()
    vi.x1 <- vi.x2 <- vi.x3 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
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
      SNCT <- SimNetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                    output = "estimation", # detection, estimation
                                    A = adiac_matrix,
                                    p = rep(probT,n), Ne = NeighNum,
                                    W = w, Y = y, X = X, M = M, G = g,
                                    mdisc = 5, mest = 5,
                                    minpopfrac = 1,
                                    depth = 2,
                                    fracpredictors = 1,
                                    minsize = 20,
                                    n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" | # change here
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.x3[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.4") |
                                               str_detect(SNCT$FILTER, "X.5") |
                                               str_detect(SNCT$FILTER, "X.6") |
                                               str_detect(SNCT$FILTER, "X.7") |
                                               str_detect(SNCT$FILTER, "X.8") |
                                               str_detect(SNCT$FILTER, "X.9") |
                                               str_detect(SNCT$FILTER, "X.10")])
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.est.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.est.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.est.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.est.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      ## Extract Values for the testimation Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                        rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                  rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1000_2 <- 0
          effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          effects.test1101_2 <- 0 
          se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1000_2 <- 0
          se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          se.test1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1000_1 <- 0
          effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.test1101_1 <- 0
          effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1000_1 <- 0
          se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1101_1 <- 0
          se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1000_1 <- 0
        effects.test1000_2 <- 0
        effects.test1101_1 <- 0
        effects.test1101_2 <- 0
        se.test1000_1 <- 0
        se.test1000_2 <- 0
        se.test1101_1 <- 0
        se.test1101_2 <- 0
      }
      
      tau.test.1000[j*2-1, which(seq==i)] <- effects.test1000_1
      tau.test.1000[j*2, which(seq==i)] <- effects.test1000_2
      se.tau.test.1000[j*2-1, which(seq==i)] <-  se.test1000_1
      se.tau.test.1000[j*2, which(seq==i)] <-  se.test1000_2
      tau.test.1101[j*2-1, which(seq==i)] <- effects.test1101_1
      tau.test.1101[j*2, which(seq==i)] <- effects.test1101_2
      se.tau.test.1101[j*2-1, which(seq==i)] <-  se.test1101_1
      se.tau.test.1101[j*2, which(seq==i)] <-  se.test1101_2
      
      ## Clear Memory
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    
    ## Return the values
    
    list(correct.rules, tau.est.1000, tau.est.1101, se.tau.est.1000, se.tau.est.1101,
         tau.test.1000, tau.test.1101, se.tau.test.1000, se.tau.test.1101, vi.x1, vi.x2, vi.3, vi.X)
  }
})

## Extract the Results 

correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
tau_est_1101 <- na.omit(matrix[[3]])
se_tau_est_1000 <- na.omit(matrix[[4]])
se_tau_est_1101 <- na.omit(matrix[[5]])
tau_test_1000 <- na.omit(matrix[[6]])
tau_test_1101 <- na.omit(matrix[[7]])
se_tau_test_1000 <- na.omit(matrix[[8]])
se_tau_test_1101 <- na.omit(matrix[[9]])
vi_x1 <- na.omit(matrix[[10]])
vi_x2 <- na.omit(matrix[[11]])
vi_x3 <- na.omit(matrix[[12]])
vi_X <- na.omit(matrix[[13]])


## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_x3 <- colMeans(vi_x3)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered

tau_est_1000[tau_est_1000==0] <- NA
tau_est_1101[tau_est_1101==0] <- NA
tau_test_1000[tau_est_1000==0] <- NA
tau_test_1101[tau_est_1101==0] <- NA


## Mean Squared Error (Estimation/Test Set)

mse_tau_est_1000 <- c()
mse_tau_test_1000 <- c()
mse_tau_est_1101 <- c()
mse_tau_test_1101 <- c()

for(i in seq){
  mse_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_tau_est_1000 <- c()
bias_tau_test_1000 <- c()
bias_tau_est_1101 <- c()
bias_tau_test_1101 <- c()

for(i in seq){
  bias_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1000 <- data.frame()
coverage_test_1000 <- data.frame()
coverage_est_1101 <- data.frame()
coverage_test_1101 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_est_1000)) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000[j,which(seq==i)] <- between(i, abs(tau_est_1000[j,which(seq==i)]) - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                    abs(tau_est_1000[j,which(seq==i)]) + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1000)) {
    if (is.na(tau_test_1000[j,which(seq==i)])==FALSE) {
      coverage_test_1000[j,which(seq==i)] <- between(i, abs(tau_test_1000[j,which(seq==i)]) - 1.96*se_tau_test_1000[j,which(seq==i)],
                                                     abs(tau_test_1000[j,which(seq==i)]) + 1.96*se_tau_test_1000[j,which(seq==i)])
    }
    else {
      coverage_test_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1000 <- colMeans(coverage_test_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_est_1101)) {
    if (is.na(tau_est_1101[j,which(seq==i)])==FALSE) {
      coverage_est_1101[j,which(seq==i)] <- between(i, abs(tau_est_1101[j,which(seq==i)]) - 1.96*se_tau_est_1101[j,which(seq==i)],
                                                    abs(tau_est_1101[j,which(seq==i)]) + 1.96*se_tau_est_1101[j,which(seq==i)])
    }
    else {
      coverage_est_1101[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1101 <- colMeans(coverage_est_1101, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1101)) {
    if (is.na(tau_test_1101[j,which(seq==i)])==FALSE) {
      coverage_test_1101[j,which(seq==i)] <- between(i, abs(tau_test_1101[j,which(seq==i)]) - 1.96*se_tau_test_1101[j,which(seq==i)],
                                                     abs(tau_test_1101[j,which(seq==i)]) + 1.96*se_tau_test_1101[j,which(seq==i)])
    }
    else {
      coverage_test_1101[j,which(seq==i)] <- NA
    }
  }  
}

## Create a Matrix for the Results

nctree <- cbind(avg_correct_rules, mse_tau_est_1000, mse_tau_est_1101, bias_tau_est_1000, bias_tau_est_1101, coverage_est_1000, coverage_est_1101,
                mse_tau_test_1000, mse_tau_test_1101, bias_tau_test_1000, bias_tau_test_1101, coverage_test_1000, coverage_test_1101, vi_x1, vi_x2, vi_3, vi_X)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_est_1000", "mse_tau_est_1101",
                      "bias_tau_est_1000", "bias_tau_est_1101",
                      "coverage_est_1000", "coverage_est_1101",
                      "mse_tau_test_1000", "mse_tau_test_1101",
                      "bias_tau_test_1000", "bias_tau_test_1101",
                      "coverage_test_1000", "coverage_test_1101",
                      "vi_x1", "vi_x2", "vi_x3", "vi_X")
write.csv(nctree, file = "two_main_effects_overlap.csv")

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
    tau.est.1000 <- data.frame()
    tau.est.1101 <- data.frame()
    se.tau.est.1000 <- data.frame()
    se.tau.est.1101 <- data.frame()
    tau.test.1000 <- data.frame()
    tau.test.1101 <- data.frame()
    se.tau.test.1000 <- data.frame()
    se.tau.test.1101 <- data.frame()
    vi.x1 <- vi.x2 <- xi.x3 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
      # Compute number of treated neighbors and consequently G_i
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
      tau1101[x1==0 & x2==0] <- i*2 # double effect size
      tau1101[x1==1 & x3==1] <- -i*2 # different HDVs & double effect size
      
      ## Generate Treatment Effects
      y01 <- rnorm(n)
      y11 <- y01 + tau1101
      y00 <- rnorm(n)
      y10 <- y00 + tau1000
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      SNCT <- SimNetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                    output = "estimation", # detection, estimation
                                    A = adiac_matrix,
                                    p = rep(probT,n), Ne = NeighNum,
                                    W = w, Y = y, X = X, M = M, G = g,
                                    mdisc = 5, mest = 5,
                                    minpopfrac = 1,
                                    depth = 2,
                                    fracpredictors = 1,
                                    minsize = 20,
                                    n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" | # change here
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.x3[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.4") |
                                               str_detect(SNCT$FILTER, "X.5") |
                                               str_detect(SNCT$FILTER, "X.6") |
                                               str_detect(SNCT$FILTER, "X.7") |
                                               str_detect(SNCT$FILTER, "X.8") |
                                               str_detect(SNCT$FILTER, "X.9") |
                                               str_detect(SNCT$FILTER, "X.10")])
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.est.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.est.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.est.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.est.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      ## Extract Values for the testimation Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                        rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | # change here
                                                  rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
        se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1000_2 <- 0
          effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          effects.test1101_2 <- 0 
          se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1000_2 <- 0
          se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] # change here
          se.test1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1000_1 <- 0
          effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.test1101_1 <- 0
          effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1000_1 <- 0
          se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1101_1 <- 0
          se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1000_1 <- 0
        effects.test1000_2 <- 0
        effects.test1101_1 <- 0
        effects.test1101_2 <- 0
        se.test1000_1 <- 0
        se.test1000_2 <- 0
        se.test1101_1 <- 0
        se.test1101_2 <- 0
      }
      
      tau.test.1000[j*2-1, which(seq==i)] <- effects.test1000_1
      tau.test.1000[j*2, which(seq==i)] <- effects.test1000_2
      se.tau.test.1000[j*2-1, which(seq==i)] <-  se.test1000_1
      se.tau.test.1000[j*2, which(seq==i)] <-  se.test1000_2
      tau.test.1101[j*2-1, which(seq==i)] <- effects.test1101_1
      tau.test.1101[j*2, which(seq==i)] <- effects.test1101_2
      se.tau.test.1101[j*2-1, which(seq==i)] <-  se.test1101_1
      se.tau.test.1101[j*2, which(seq==i)] <-  se.test1101_2
      
      ## Clear Memory
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    
    ## Return the values
    
    list(correct.rules, tau.est.1000, tau.est.1101, se.tau.est.1000, se.tau.est.1101,
         tau.test.1000, tau.test.1101, se.tau.test.1000, se.tau.test.1101, vi.x1, vi.x2, vi.x3, vi.X)
  }
})

## Extract the Results 

correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
tau_est_1101 <- na.omit(matrix[[3]])
se_tau_est_1000 <- na.omit(matrix[[4]])
se_tau_est_1101 <- na.omit(matrix[[5]])
tau_test_1000 <- na.omit(matrix[[6]])
tau_test_1101 <- na.omit(matrix[[7]])
se_tau_test_1000 <- na.omit(matrix[[8]])
se_tau_test_1101 <- na.omit(matrix[[9]])
vi_x1 <- na.omit(matrix[[10]])
vi_x2 <- na.omit(matrix[[11]])
vi_x3 <- na.omit(matrix[[12]])
vi_X <- na.omit(matrix[[13]])

## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_x3 <- colMeans(vi_x3)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered

tau_est_1000[tau_est_1000==0] <- NA
tau_est_1101[tau_est_1101==0] <- NA
tau_test_1000[tau_est_1000==0] <- NA
tau_test_1101[tau_est_1101==0] <- NA


## Mean Squared Error (Estimation/Test Set)

mse_tau_est_1000 <- c()
mse_tau_test_1000 <- c()
mse_tau_est_1101 <- c()
mse_tau_test_1101 <- c()

for(i in seq){
  mse_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_tau_est_1000 <- c()
bias_tau_test_1000 <- c()
bias_tau_est_1101 <- c()
bias_tau_test_1101 <- c()

for(i in seq){
  bias_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1000 <- data.frame()
coverage_test_1000 <- data.frame()
coverage_est_1101 <- data.frame()
coverage_test_1101 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_est_1000)) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000[j,which(seq==i)] <- between(i, abs(tau_est_1000[j,which(seq==i)]) - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                    abs(tau_est_1000[j,which(seq==i)]) + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1000)) {
    if (is.na(tau_test_1000[j,which(seq==i)])==FALSE) {
      coverage_test_1000[j,which(seq==i)] <- between(i, abs(tau_test_1000[j,which(seq==i)]) - 1.96*se_tau_test_1000[j,which(seq==i)],
                                                     abs(tau_test_1000[j,which(seq==i)]) + 1.96*se_tau_test_1000[j,which(seq==i)])
    }
    else {
      coverage_test_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1000 <- colMeans(coverage_test_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_est_1101)) {
    if (is.na(tau_est_1101[j,which(seq==i)])==FALSE) {
      coverage_est_1101[j,which(seq==i)] <- between(i, abs(tau_est_1101[j,which(seq==i)]) - 1.96*se_tau_est_1101[j,which(seq==i)],
                                                    abs(tau_est_1101[j,which(seq==i)]) + 1.96*se_tau_est_1101[j,which(seq==i)])
    }
    else {
      coverage_est_1101[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1101 <- colMeans(coverage_est_1101, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1101)) {
    if (is.na(tau_test_1101[j,which(seq==i)])==FALSE) {
      coverage_test_1101[j,which(seq==i)] <- between(i, abs(tau_test_1101[j,which(seq==i)]) - 1.96*se_tau_test_1101[j,which(seq==i)],
                                                     abs(tau_test_1101[j,which(seq==i)]) + 1.96*se_tau_test_1101[j,which(seq==i)])
    }
    else {
      coverage_test_1101[j,which(seq==i)] <- NA
    }
  }  
}

## Create a Matrix for the Results

nctree <- cbind(avg_correct_rules, mse_tau_est_1000, mse_tau_est_1101, bias_tau_est_1000, bias_tau_est_1101, coverage_est_1000, coverage_est_1101,
                mse_tau_test_1000, mse_tau_test_1101, bias_tau_test_1000, bias_tau_test_1101, coverage_test_1000, coverage_test_1101, vi_x1, vi_x2, vi_x3, vi.X)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_est_1000", "mse_tau_est_1101",
                      "bias_tau_est_1000", "bias_tau_est_1101",
                      "coverage_est_1000", "coverage_est_1101",
                      "mse_tau_test_1000", "mse_tau_test_1101",
                      "bias_tau_test_1000", "bias_tau_test_1101",
                      "coverage_test_1000", "coverage_test_1101",
                      "vi_x1", "vi_x2", "vi_x3", "vi_X")
write.csv(nctree, file = "two_main_effects_overlap_effect_size.csv")


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
    tau.est.1000 <- data.frame()
    tau.est.1101 <- data.frame()
    se.tau.est.1000 <- data.frame()
    se.tau.est.1101 <- data.frame()
    tau.test.1000 <- data.frame()
    tau.test.1101 <- data.frame()
    se.tau.test.1000 <- data.frame()
    se.tau.test.1101 <- data.frame()
    vi.x1 <- vi.x2 <- vi.X <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency matrix
      adiac_matrix <- genmultnet(N=N, m=15, method="er", param = 0.05)
      
      # Group Indicator 
      M <- c(rep(1:m, N/m))
      M <- sort(M)
      levels(M) <- c(1:m)
      
      # Generate treatment
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
      
      # Degree
      neigh <- rowSums(adiac_matrix)
      
      # Compute number of treated neighbors and consequently G_i
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
      SNCT <- SimNetworkCausalTrees(effweights = c(0.5,0.5,0,0), method = "composite", # singular
                                    output = "estimation", # detection, estimation
                                    A = adiac_matrix,
                                    p = rep(probT,n), Ne = NeighNum,
                                    W = w, Y = y, X = X, M = M, G = g,
                                    mdisc = 5, mest = 5,
                                    minpopfrac = 1,
                                    depth = 2,
                                    fracpredictors = 1,
                                    minsize = 20,
                                    n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      
      ## Extract the Correct Rules
      
      correct <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Variables Importance
      
      vi.x1[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.1")])
      vi.x2[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.2")])
      vi.X[j, which(seq==i)] <- sum(SNCT$GOF[str_detect(SNCT$FILTER, "X.3") | 
                                               str_detect(SNCT$FILTER, "X.4") |
                                               str_detect(SNCT$FILTER, "X.5") |
                                               str_detect(SNCT$FILTER, "X.6") |
                                               str_detect(SNCT$FILTER, "X.7") |
                                               str_detect(SNCT$FILTER, "X.8") |
                                               str_detect(SNCT$FILTER, "X.9") |
                                               str_detect(SNCT$FILTER, "X.10")])
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      
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
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      tau.est.1101[j*2-1, which(seq==i)] <- effects.est1101_1
      tau.est.1101[j*2, which(seq==i)] <- effects.est1101_2
      se.tau.est.1101[j*2-1, which(seq==i)] <-  se.est1101_1
      se.tau.est.1101[j*2, which(seq==i)] <-  se.est1101_2
      
      ## Extract Values for the testimation Set (Get 0 if the Causal Rule was not identified)
      
      if (correct==2) {
        effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                        rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                        rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
        se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" |
                                                  rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
        se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.2<1 & data_tree$X.1<1" |
                                                  rule.sel=="data_tree$X.1<1 & data_tree$X.2<1")]
      }
      
      if (correct==1){
        
        if (length((which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")))==1){
          effects.test1000_1 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1000_2 <- 0
          effects.test1101_1 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          effects.test1101_2 <- 0
          se.test1000_1 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1000_2 <- 0
          se.test1101_1 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1>=1")]
          se.test1101_2 <- 0
        }
        
        if (length((which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")))==1){
          effects.test1000_1 <- 0
          effects.test1000_2 <- SNCT$EFF1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          effects.test1101_1 <- 0
          effects.test1101_2 <- SNCT$EFF1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1000_1 <- 0
          se.test1000_2 <- SNCT$SE1000_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
          se.test1101_1 <- 0
          se.test1101_2 <- SNCT$SE1101_TEST[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")]
        }
        
      }  
      
      if (correct==0) {
        effects.test1000_1 <- 0
        effects.test1000_2 <- 0
        effects.test1101_1 <- 0
        effects.test1101_2 <- 0
        se.test1000_1 <- 0
        se.test1000_2 <- 0
        se.test1101_1 <- 0
        se.test1101_2 <- 0
      }
      
      tau.test.1000[j*2-1, which(seq==i)] <- effects.test1000_1
      tau.test.1000[j*2, which(seq==i)] <- effects.test1000_2
      se.tau.test.1000[j*2-1, which(seq==i)] <-  se.test1000_1
      se.tau.test.1000[j*2, which(seq==i)] <-  se.test1000_2
      tau.test.1101[j*2-1, which(seq==i)] <- effects.test1101_1
      tau.test.1101[j*2, which(seq==i)] <- effects.test1101_2
      se.tau.test.1101[j*2-1, which(seq==i)] <-  se.test1101_1
      se.tau.test.1101[j*2, which(seq==i)] <-  se.test1101_2
      
      ## Clear Memory
      
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est1101_1, effects.est1101_2,
         se.est1000_1, se.est1000_2,
         se.est1101_1, se.est1101_2)
    }
    
    ## Return the values
    
    list(correct.rules, tau.est.1000, tau.est.1101, se.tau.est.1000, se.tau.est.1101,
         tau.test.1000, tau.test.1101, se.tau.test.1000, se.tau.test.1101, vi.x1, vi.x2, vi.X)
  }
})

## Extract the Results 

correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
tau_est_1101 <- na.omit(matrix[[3]])
se_tau_est_1000 <- na.omit(matrix[[4]])
se_tau_est_1101 <- na.omit(matrix[[5]])
tau_test_1000 <- na.omit(matrix[[6]])
tau_test_1101 <- na.omit(matrix[[7]])
se_tau_test_1000 <- na.omit(matrix[[8]])
se_tau_test_1101 <- na.omit(matrix[[9]])
vi_x1 <- na.omit(matrix[[10]])
vi_x2 <- na.omit(matrix[[11]])
vi_X <- na.omit(matrix[[12]])

## Correct Rules

avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 

vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered

tau_est_1000[tau_est_1000==0] <- NA
tau_est_1101[tau_est_1101==0] <- NA
tau_test_1000[tau_est_1000==0] <- NA
tau_test_1101[tau_est_1101==0] <- NA


## Mean Squared Error (Estimation/Test Set)

mse_tau_est_1000 <- c()
mse_tau_test_1000 <- c()
mse_tau_est_1101 <- c()
mse_tau_test_1101 <- c()

for(i in seq){
  mse_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

for(i in seq){
  mse_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i)^2 , na.rm = TRUE )
}

## Bias (Estimation/Test Set)

bias_tau_est_1000 <- c()
bias_tau_test_1000 <- c()
bias_tau_est_1101 <- c()
bias_tau_test_1101 <- c()

for(i in seq){
  bias_tau_est_1000[which(seq==i)] <- mean( (abs(tau_est_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1000[which(seq==i)] <- mean( (abs(tau_test_1000[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_est_1101[which(seq==i)] <- mean( (abs(tau_est_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

for(i in seq){
  bias_tau_test_1101[which(seq==i)] <- mean( (abs(tau_test_1101[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage (Estimation/Test Set)

coverage_est_1000 <- data.frame()
coverage_test_1000 <- data.frame()
coverage_est_1101 <- data.frame()
coverage_test_1101 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_est_1000)) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000[j,which(seq==i)] <- between(i, abs(tau_est_1000[j,which(seq==i)]) - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                    abs(tau_est_1000[j,which(seq==i)]) + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1000)) {
    if (is.na(tau_test_1000[j,which(seq==i)])==FALSE) {
      coverage_test_1000[j,which(seq==i)] <- between(i, abs(tau_test_1000[j,which(seq==i)]) - 1.96*se_tau_test_1000[j,which(seq==i)],
                                                     abs(tau_test_1000[j,which(seq==i)]) + 1.96*se_tau_test_1000[j,which(seq==i)])
    }
    else {
      coverage_test_1000[j,which(seq==i)] <- NA
    }
  }  
}
coverage_test_1000 <- colMeans(coverage_test_1000, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_est_1101)) {
    if (is.na(tau_est_1101[j,which(seq==i)])==FALSE) {
      coverage_est_1101[j,which(seq==i)] <- between(i, abs(tau_est_1101[j,which(seq==i)]) - 1.96*se_tau_est_1101[j,which(seq==i)],
                                                    abs(tau_est_1101[j,which(seq==i)]) + 1.96*se_tau_est_1101[j,which(seq==i)])
    }
    else {
      coverage_est_1101[j,which(seq==i)] <- NA
    }
  }  
}
coverage_est_1101 <- colMeans(coverage_est_1101, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_test_1101)) {
    if (is.na(tau_test_1101[j,which(seq==i)])==FALSE) {
      coverage_test_1101[j,which(seq==i)] <- between(i, abs(tau_test_1101[j,which(seq==i)]) - 1.96*se_tau_test_1101[j,which(seq==i)],
                                                     abs(tau_test_1101[j,which(seq==i)]) + 1.96*se_tau_test_1101[j,which(seq==i)])
    }
    else {
      coverage_test_1101[j,which(seq==i)] <- NA
    }
  }  
}

## Create a Matrix for the Results

nctree <- cbind(avg_correct_rules, mse_tau_est_1000, mse_tau_est_1101, bias_tau_est_1000, bias_tau_est_1101, coverage_est_1000, coverage_est_1101,
                mse_tau_test_1000, mse_tau_test_1101, bias_tau_test_1000, bias_tau_test_1101, coverage_test_1000, coverage_test_1101, vi.x1, vi.x2, vi.X)
colnames(nctree) <- c("correct_rules",
                      "mse_tau_est_1000", "mse_tau_est_1101",
                      "bias_tau_est_1000", "bias_tau_est_1101",
                      "coverage_est_1000", "coverage_est_1101",
                      "mse_tau_test_1000", "mse_tau_test_1101",
                      "bias_tau_test_1000", "bias_tau_test_1101",
                      "coverage_test_1000", "coverage_test_1101",
                      "vi_x1", "vi_x2", "vi_X")
write.csv(nctree, file = "two_main_effects_correlated.csv")

stopCluster()

## End of Simulations
