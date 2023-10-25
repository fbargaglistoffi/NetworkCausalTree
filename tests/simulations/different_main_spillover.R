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
library(dplyr)

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
# gsize: size of each cluster
# param: link probability
N = 3000
m = 30
n_cov = 10
prob = 0.5 
gsize = N/m
param = 0.01
mu = rep(0, n_cov)
rho = 0
seq <- seq(0.1, 10.1, 1)
nsim = 500

# Set up Parallel Computation
# Setup parallel backend to use many processors
cl <- makeCluster(20)
registerDoParallel(cl)

# This function takes an arbitrary number of lists x all of which much have the same structure    
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}


##################################
##    4 RULES & MAIN + SPILL    ##
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
    correct.rules <-  data.frame()
    tau.est.1000 <- data.frame()
    eta.est.0100 <- data.frame()
    se.tau.est.1000 <- data.frame()
    se.eta.est.0100 <- data.frame()
    vi.x1 <- vi.x2 <- vi.x3 <- vi.X <- data.frame()
    correct.rules.main <- data.frame()
    correct.rules.spil <- data.frame()
    correct.rules.composite.main <- correct.rules.composite.spil <- data.frame()
    correct.singlemain.spil <- correct.singlespil.main <- data.frame()
    
    for (i in seq)   {
      
      ###########################
      # Data Generating Process #
      ##########################
      
      # Adjacency Matrix
      adiac_matrix <- genmultnet(N=N, m=m, method="er", param = param)
      
      # Group Indicator 
      M=c(rep(1:m,gsize))
      M=sort(M)
      levels(M)<-c(1:m)
      
      # Generate treatment (keep prob fixed at 0.4)
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
      
      # Neighbors
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
      tau1000[x1==0 & x2==0] <- 3*i
      tau1000[x1==0 & x2==1] <- i
      
      # Tau 0100
      eta0100 <- rep(0, n)
      eta0100[x1==1 & x3==0] <- 3*i
      eta0100[x1==1 & x3==1] <- i
      
      ## Generate Treatment Effects
      y00 <- rnorm(n)
      y11 <- rnorm(n)
      y10 <- y00 + tau1000
      y01 <- y00 + eta0100
      
      ## Generate Outcome
      y <- y00*(1-w)*(1-g) + y10*w*(1-g) + y01*(1-w)*g + y11*w*g
      
      
      ## Run the  whole function
      
      ## Composite Tree
      SNCT <- NetworkCausalTrees(effweights = c(0.5,0,0,0.5), method = "composite", # singular
                                 output = "estimation", # detection, estimation
                                 A = adiac_matrix,
                                 p = rep(probT,n), Ne = NeighNum,
                                 W = w, Y = y, X = X, M = M, G = g,
                                 mdisc = m/2, mest = m/2,
                                 minpopfrac = 1,
                                 depth = 2,
                                 fracpredictors = 1,
                                 minsize= N/100,
                                 n_trees = 1) 
      
      ## Tree for Main Treatment Effect
      SNCT_main <- NetworkCausalTrees(effweights = c(1,0,0,0), method = "singular",
                                      output = "estimation", # detection, estimation
                                      A = adiac_matrix,
                                      p = rep(probT,n), Ne = NeighNum,
                                      W = w, Y = y, X = X, M = M, G = g,
                                      mdisc = m/2, mest = m/2,
                                      minpopfrac = 1,
                                      depth = 2,
                                      fracpredictors = 1,
                                      minsize= N/100,
                                      n_trees = 1) 
      
      ## Tree for Spillover Effect
      SNCT_spil <- NetworkCausalTrees(effweights = c(0,0,0,1), method = "singular",
                                      output = "estimation", # detection, estimation
                                      A = adiac_matrix,
                                      p = rep(probT,n), Ne = NeighNum,
                                      W = w, Y = y, X = X, M = M, G = g,
                                      mdisc = m/2, mest = m/2,
                                      minpopfrac = 1,
                                      depth = 2,
                                      fracpredictors = 1,
                                      minsize= N/100,
                                      n_trees = 1) 
      
      rule.sel <- SNCT$FILTER
      rule.main <- SNCT_main$FILTER
      rule.spil <- SNCT_spil$FILTER
      
      ############################
      ## Composite Tree results ##
      ############################
      
      ## Extract the Correct Rules
      correct <- length(which(rule.sel=="data_tree$X.1<1 & data_tree$X.2>=1" |
                                rule.sel=="data_tree$X.2>=1 & data_tree$X.1<1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | 
                                rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" |
                                rule.sel=="data_tree$X.1>=1 & data_tree$X.3<1" | 
                                rule.sel=="data_tree$X.3<1 & data_tree$X.1>=1" | 
                                rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules[j, which(seq==i)] <- correct
      
      ## Number rules
      SNCT$RULE.NUM <- rep(0, nrow(SNCT))
      SNCT$RULE.NUM[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" | rule.sel=="data_tree$X.2<1 & data_tree$X.1<1")] <- 1
      SNCT$RULE.NUM[which(rule.sel=="data_tree$X.1<1 & data_tree$X.2>=1" | rule.sel=="data_tree$X.2>=1 & data_tree$X.1<1")] <- 2
      SNCT$RULE.NUM[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1")] <-3
      SNCT$RULE.NUM[which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3<1" | rule.sel=="data_tree$X.3<1 & data_tree$X.1>=1")] <-4
      
      ## Variables Importance
      vi.x1[j, which(seq==i)] <- str_detect(SNCT$FILTER, "X.1")[2] # the second element corresponds to the first split
      vi.x2[j, which(seq==i)] <- str_detect(SNCT$FILTER, "X.2")[2]
      vi.x3[j, which(seq==i)] <- str_detect(SNCT$FILTER, "X.3")[2]
      vi.X[j, which(seq==i)] <-  cbind(str_detect(SNCT$FILTER, "X.4")[2] |
                                         str_detect(SNCT$FILTER, "X.5")[2] |  
                                         str_detect(SNCT$FILTER, "X.6")[2] |
                                         str_detect(SNCT$FILTER, "X.7")[2] |  
                                         str_detect(SNCT$FILTER, "X.8")[2] |
                                         str_detect(SNCT$FILTER, "X.9")[2] |  
                                         str_detect(SNCT$FILTER, "X.10")[2])
      
      ## Extract Values for the Estimation Set (Get 0 if the Causal Rule was not identified)
      if (correct==4) {
        effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)] 
        effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
        effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)] 
        effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
        se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)] 
        se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
        se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)] 
        se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
      }
      
      if (correct==3) {
        
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==2 | SNCT$RULE.NUM==3)))==3){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)] 
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- 0 
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)] 
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- 0
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==2 | SNCT$RULE.NUM==4)))==3){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)] 
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)] 
          effects.est0100_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)] 
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)] 
          se.est0100_2 <- 0
        }
        
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==3 | SNCT$RULE.NUM==4)))==3){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)] 
          effects.est1000_2 <- 0
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)] 
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)] 
          se.est1000_2 <- 0
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)] 
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
        if (length((which(SNCT$RULE.NUM==2 | SNCT$RULE.NUM==3 | SNCT$RULE.NUM==4)))==3){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)] 
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)] 
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
      }
      
      if (correct==2){ 
        # tau1000_1 & tau1000_2
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==2)))==2){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)]
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)] 
          effects.est0100_1 <- 0
          effects.est0100_2 <- 0 
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)]
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)] 
          se.est0100_1 <- 0
          se.est0100_2 <- 0
        }
        
        # tau1000_1 & tau0100_2
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==3)))==2){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)]
          effects.est1000_2 <- 0
          effects.est0100_1 <- 0
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)] 
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)]
          se.est1000_2 <- 0
          se.est0100_1 <- 0 
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
        # tau1000_2 & tau0100_2
        if (length((which(SNCT$RULE.NUM==2 | SNCT$RULE.NUM==3)))==2){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- 0
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- 0
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
        # tau1000_1 & tau0100_1
        if (length((which(SNCT$RULE.NUM==1 | SNCT$RULE.NUM==4)))==2){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)]
          effects.est1000_2 <- 0
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)]
          effects.est0100_2 <- 0
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)]
          se.est1000_2 <- 0
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)]
          se.est0100_2 <- 0
        }
        
        # tau1000_2 & tau0100_1
        if (length((which(SNCT$RULE.NUM==2 | SNCT$RULE.NUM==4)))==2){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)]
          effects.est0100_2 <- 0
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)]
          se.est0100_2 <- 0
        }
        
        # tau0100_1 & tau0100_2
        if (length((which(SNCT$RULE.NUM==3 | SNCT$RULE.NUM==4)))==2){
          effects.est1000_1 <- 0
          effects.est1000_2 <- 0
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)]
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- 0
          se.est1000_2 <- 0
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)]
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
      }  
      
      if (correct==1){ 
        # tau1000_1 
        if (length((which(SNCT$RULE.NUM==1)))==1){
          effects.est1000_1 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==1)]
          effects.est1000_2 <- 0
          effects.est0100_1 <- 0
          effects.est0100_2 <- 0 
          se.est1000_1 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==1)]
          se.est1000_2 <- 0
          se.est0100_1 <- 0
          se.est0100_2 <- 0
        }
        
        # tau1000_2
        if (length((which(SNCT$RULE.NUM==2)))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- SNCT$EFF1000_EST[which(SNCT$RULE.NUM==2)]
          effects.est0100_1 <- 0
          effects.est0100_2 <- 0 
          se.est1000_1 <- 0
          se.est1000_2 <- SNCT$SE1000_EST[which(SNCT$RULE.NUM==2)]
          se.est0100_1 <- 0
          se.est0100_2 <- 0
        }
        
        # tau0100_2
        if (length((which(SNCT$RULE.NUM==3)))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- 0
          effects.est0100_1 <- 0
          effects.est0100_2 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==3)]
          se.est1000_1 <- 0
          se.est1000_2 <- 0
          se.est0100_1 <- 0
          se.est0100_2 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==3)]
        }
        
        # tau0100_1
        if (length((which(SNCT$RULE.NUM==4)))==1){
          effects.est1000_1 <- 0
          effects.est1000_2 <- 0
          effects.est0100_1 <- SNCT$EFF0100_EST[which(SNCT$RULE.NUM==4)]
          effects.est0100_2 <- 0
          se.est1000_1 <- 0
          se.est1000_2 <- 0
          se.est0100_1 <- SNCT$SE0100_EST[which(SNCT$RULE.NUM==4)]
          se.est0100_2 <- 0
        }
      }  
      
      if (correct==0) {
        effects.est1000_1 <- 0
        effects.est1000_2 <- 0
        effects.est0100_1 <- 0
        effects.est0100_2 <- 0
        se.est1000_1 <- 0
        se.est1000_2 <- 0
        se.est0100_1 <- 0
        se.est0100_2 <- 0
      }
      
      tau.est.1000[j*2-1, which(seq==i)] <- effects.est1000_1
      tau.est.1000[j*2, which(seq==i)] <- effects.est1000_2
      se.tau.est.1000[j*2-1, which(seq==i)] <-  se.est1000_1
      se.tau.est.1000[j*2, which(seq==i)] <-  se.est1000_2
      eta.est.0100[j*2-1, which(seq==i)] <- effects.est0100_1
      eta.est.0100[j*2, which(seq==i)] <- effects.est0100_2
      se.eta.est.0100[j*2-1, which(seq==i)] <-  se.est0100_1
      se.eta.est.0100[j*2, which(seq==i)] <-  se.est0100_2
      
      #############################
      ##  Treatment Tree Results ## 
      #############################
      
      ## Extract the Correct Rules (Composite Tree)
      correct.composite.main <- length(which(rule.sel=="data_tree$X.1<1 & data_tree$X.2>=1" |
                                               rule.sel=="data_tree$X.2>=1 & data_tree$X.1<1" |
                                               rule.sel=="data_tree$X.1<1 & data_tree$X.2<1" |
                                               rule.sel=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules.composite.main[j, which(seq==i)] <- correct.composite.main
      
      ## Extract the Correct Rules (Spillover Tree)
      correct.singlespil <- length(which(rule.spil=="data_tree$X.1<1 & data_tree$X.2>=1" |
                                           rule.spil=="data_tree$X.2>=1 & data_tree$X.1<1" |
                                           rule.spil=="data_tree$X.1<1 & data_tree$X.2<1" |
                                           rule.spil=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.singlespil.main[j, which(seq==i)] <- correct.singlespil
      
      ## Extract the Correct Rules (Treatment Tree)
      correct.main <- length(which(rule.main=="data_tree$X.1<1 & data_tree$X.2>=1" |
                                     rule.main=="data_tree$X.2>=1 & data_tree$X.1<1" |
                                     rule.main=="data_tree$X.1<1 & data_tree$X.2<1" |
                                     rule.main=="data_tree$X.2<1 & data_tree$X.1<1"))
      correct.rules.main[j, which(seq==i)] <- correct.main
      
      
      #############################
      ##  Spillover Tree Results  # 
      #############################
      
      ## Extract the Correct Rules (Composite Tree)
      correct.composite.spil <- length(which(rule.sel=="data_tree$X.1>=1 & data_tree$X.3>=1" | 
                                               rule.sel=="data_tree$X.3>=1 & data_tree$X.1>=1" | 
                                               rule.sel=="data_tree$X.1>=1 & data_tree$X.3<1" |
                                               rule.sel=="data_tree$X.3<1 & data_tree$X.1>=1"))
      correct.rules.composite.spil[j, which(seq==i)] <- correct.composite.spil
      
      ## Extract the Correct Rules (Treatment Tree)
      correct.singlemain <- length(which(rule.main=="data_tree$X.1>=1 & data_tree$X.3>=1" | 
                                           rule.main=="data_tree$X.3>=1 & data_tree$X.1>=1" | 
                                           rule.main=="data_tree$X.1<1 & data_tree$X.3<1" |
                                           rule.main=="data_tree$X.3<1 & data_tree$X.1<1"))
      correct.singlemain.spil[j, which(seq==i)] <- correct.singlemain
      
      ## Extract the Correct Rules (Spillover Tree)
      correct.spil <- length(which(rule.spil=="data_tree$X.1>=1 & data_tree$X.3>=1" | 
                                     rule.spil=="data_tree$X.3>=1 & data_tree$X.1>=1" | 
                                     rule.spil=="data_tree$X.1>=1 & data_tree$X.3<1" |
                                     rule.spil=="data_tree$X.3<1 & data_tree$X.1>=1"))
      correct.rules.spil[j, which(seq==i)] <- correct.spil
      
      ## Clear Memory
      rm(w, g, y, M, X, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
         probT, NeighNum, n, adiac_matrix, rule.sel, SNCT,
         effects.est1000_1, effects.est1000_2,
         effects.est0100_1, effects.est0100_2,
         se.est1000_1, se.est1000_2,
         se.est0100_1, se.est0100_2)
    }
    
    ## Return the values
    list(correct.rules, tau.est.1000, eta.est.0100, se.tau.est.1000, se.eta.est.0100, vi.x1, vi.x2, vi.x3, vi.X,
         correct.rules.main, 
         correct.rules.spil, 
         correct.rules.composite.main,
         correct.rules.composite.spil,
         correct.singlespil.main,
         correct.singlemain.spil)
    
  }
})

###########################
##   Store the Results   ##
###########################

###############
## Composite ##
###############

## Extract the Results 
correct_rules <- na.omit(matrix[[1]])
tau_est_1000 <- na.omit(matrix[[2]])
eta_est_0100 <- na.omit(matrix[[3]])
se_tau_est_1000 <- na.omit(matrix[[4]])
se_eta_est_0100 <- na.omit(matrix[[5]])
vi_x1 <- na.omit(matrix[[6]])
vi_x2 <- na.omit(matrix[[7]])
vi_x3 <- na.omit(matrix[[8]])
vi_X <- na.omit(matrix[[9]])

## Correct Rules
avg_correct_rules <- colMeans(correct_rules)

## Variable Importance 
vi_x1 <- colMeans(vi_x1)
vi_x2 <- colMeans(vi_x2)
vi_x3 <- colMeans(vi_x3)
vi_X <- colMeans(vi_X)

## Exclude Rules that were not discovered
tau_est_1000[tau_est_1000==0] <- NA
eta_est_0100[eta_est_0100==0] <- NA
se_tau_est_1000[se_tau_est_1000==0] <- NA 
se_eta_est_0100[se_eta_est_0100==0] <- NA 

odd_indexes<-seq(1,1000,2) 
even_indexes<-seq(2,1000,2)

## Average Effects and Average SD
avg_tau_1000_1 <-  avg_tau_1000_2 <- c()
for(i in seq){
  avg_tau_1000_1[which(seq==i)] <- mean( tau_est_1000[odd_indexes,which(seq==i)] , na.rm = TRUE )
  avg_tau_1000_2[which(seq==i)] <- mean( tau_est_1000[even_indexes,which(seq==i)] , na.rm = TRUE )
}

avg_se_1000_1 <-  avg_se_1000_2 <- c()
for(i in seq){
  avg_se_1000_1[which(seq==i)] <- mean( se_tau_est_1000[odd_indexes,which(seq==i)] , na.rm = TRUE )
  avg_se_1000_2[which(seq==i)] <- mean( se_tau_est_1000[even_indexes,which(seq==i)] , na.rm = TRUE )
}

avg_eta_0100_1 <-  avg_eta_0100_2 <- c()
for(i in seq){
  avg_eta_0100_1[which(seq==i)] <- mean( eta_est_0100[odd_indexes,which(seq==i)] , na.rm = TRUE )
  avg_eta_0100_2[which(seq==i)] <- mean( eta_est_0100[even_indexes,which(seq==i)] , na.rm = TRUE )
}

avg_se_0100_1 <-  avg_se_0100_2 <- c()
for(i in seq){
  avg_se_0100_1[which(seq==i)] <- mean( se_eta_est_0100[odd_indexes,which(seq==i)] , na.rm = TRUE )
  avg_se_0100_2[which(seq==i)] <- mean( se_eta_est_0100[even_indexes,which(seq==i)] , na.rm = TRUE )
}

## Mean Squared Error & Bias
mse_tau_est_1000 <- c()
mse_eta_est_0100 <- c()
bias_tau_est_1000 <- c()
bias_eta_est_0100 <- c()

for(i in seq){
  mse_tau_est_1000[which(seq==i)] <-  mean( c(mean( ( tau_est_1000[odd_indexes,which(seq==i)] - 3*i)^2 , na.rm = TRUE ), mean( ( tau_est_1000[even_indexes,which(seq==i)] - i)^2 , na.rm = TRUE )))
  bias_tau_est_1000[which(seq==i)] <-  mean( c(mean( ( tau_est_1000[odd_indexes,which(seq==i)] - 3*i) , na.rm = TRUE ), mean( ( tau_est_1000[even_indexes,which(seq==i)] - i) , na.rm = TRUE )))
  mse_eta_est_0100[which(seq==i)] <-  mean( c(mean( ( eta_est_0100[odd_indexes,which(seq==i)] - 3*i)^2 , na.rm = TRUE ), mean( ( eta_est_0100[even_indexes,which(seq==i)] - i)^2 , na.rm = TRUE )))
  bias_eta_est_0100[which(seq==i)] <- mean( c(mean( ( eta_est_0100[odd_indexes,which(seq==i)] - 3*i) , na.rm = TRUE ), mean( ( eta_est_0100[even_indexes,which(seq==i)] - i) , na.rm = TRUE )))
}

## Coverage
coverage_est_1000_1 <- coverage_est_1000_2 <- data.frame()

for(i in seq){
  for(j in odd_indexes) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000_1[which(odd_indexes==j),which(seq==i)] <- between(3*i, tau_est_1000[j,which(seq==i)] - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                                          tau_est_1000[j,which(seq==i)] + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000_1[which(odd_indexes==j),which(seq==i)] <- NA
    }
  }  
}

for(i in seq){
  for(j in even_indexes) {
    if (is.na(tau_est_1000[j,which(seq==i)])==FALSE) {
      coverage_est_1000_2[which(even_indexes==j),which(seq==i)] <- between(i, tau_est_1000[j,which(seq==i)] - 1.96*se_tau_est_1000[j,which(seq==i)],
                                                                           tau_est_1000[j,which(seq==i)] + 1.96*se_tau_est_1000[j,which(seq==i)])
    }
    else {
      coverage_est_1000_2[which(even_indexes==j),which(seq==i)] <- NA
    }
  }  
}

coverage_est_1000 <- rbind(coverage_est_1000_1, coverage_est_1000_2)
coverage_est_1000 <- colMeans(coverage_est_1000, na.rm = TRUE) 

coverage_est_0100_1 <- coverage_est_0100_2 <- data.frame()

for(i in seq){
  for(j in odd_indexes) {
    if (is.na(eta_est_0100[j,which(seq==i)])==FALSE) {
      coverage_est_0100_1[which(odd_indexes==j),which(seq==i)] <- between(3*i, eta_est_0100[j,which(seq==i)] - 1.96*se_eta_est_0100[j,which(seq==i)],
                                                                          eta_est_0100[j,which(seq==i)] + 1.96*se_eta_est_0100[j,which(seq==i)])
    }
    else {
      coverage_est_0100_1[which(odd_indexes==j),which(seq==i)] <- NA
    }
  }  
}

for(i in seq){
  for(j in even_indexes) {
    if (is.na(eta_est_0100[j,which(seq==i)])==FALSE) {
      coverage_est_0100_2[which(even_indexes==j),which(seq==i)] <- between(i, eta_est_0100[j,which(seq==i)] - 1.96*se_eta_est_0100[j,which(seq==i)],
                                                                           eta_est_0100[j,which(seq==i)] + 1.96*se_eta_est_0100[j,which(seq==i)])
    }
    else {
      coverage_est_0100_2[which(even_indexes==j),which(seq==i)] <- NA
    }
  }  
}

coverage_est_0100 <- rbind(coverage_est_0100_1, coverage_est_0100_2)
coverage_est_0100 <- colMeans(coverage_est_0100, na.rm = TRUE) 


## Create a Matrix for the Results
results_nctree <- cbind(avg_correct_rules,
                        mse_tau_est_1000, bias_tau_est_1000,
                        avg_tau_1000_1, avg_se_1000_1,
                        avg_tau_1000_2, avg_se_1000_2,
                        coverage_est_1000,
                        mse_eta_est_0100, bias_eta_est_0100,
                        avg_eta_0100_1, avg_se_0100_1,
                        avg_eta_0100_2, avg_se_0100_2,
                        coverage_est_0100,
                        vi_x1, vi_x2, vi_X)

colnames(results_nctree) <- c("correct_rules",
                              "MSE tau(10,00)", 
                              "Bias tau(10,00)",
                              "Rule 1 tau(10,00)",
                              "Rule 1 se(tau(10,00))", 
                              "Rule 2 tau(10,00)",
                              "Rule 2 se(tau(10,00))",
                              "Coverage tau(10,00)",
                              "MSE eta(01,00)", 
                              "Bias eta(01,00)",
                              "Rule 1 eta(01,00)",
                              "Rule 1 se(eta(01,00))",
                              "Rule 2 eta(01,00)",
                              "Rule 2 se(eta(01,00))",
                              "Coverage eta(01,00)",
                              "VI x1", "VI x2", "VI X")
write.csv(results_nctree, file = "two_main_spillover_four_effects_composite_3000.csv")

#####################
## Singular (Main) ##
#####################

## Extract Rules
correct_rules_main <- na.omit(matrix[[10]])
correct_rules_composite <- na.omit(matrix[[12]])
correct_rules_spil <- na.omit(matrix[[14]])


## Correct Rules
avg_correct_main <- colMeans(correct_rules_main)
avg_correct_composite <- colMeans(correct_rules_composite)
avg_correct_spil <- colMeans(correct_rules_spil)

results_nctree <- cbind(avg_correct_main,
                        avg_correct_composite,
                        avg_correct_spil)
colnames(results_nctree) <- c("correct_main",
                              "correct_composite",
                              "correct_spil")
write.csv(results_nctree, file = "two_main_spillover_four_effects_singular_main_3000.csv")

##########################
## Singular (Spillover) ##
##########################

## Extract Rules
correct_rules_spil <- na.omit(matrix[[11]])
correct_rules_composite <- na.omit(matrix[[13]])
correct_rules_main <- na.omit(matrix[[15]])


## Correct Rules
avg_correct_spil <- colMeans(correct_rules_spil)
avg_correct_composite <- colMeans(correct_rules_composite)
avg_correct_main <- colMeans(correct_rules_main)

results_nctree <- cbind(avg_correct_main,
                        avg_correct_composite,
                        avg_correct_spil)
colnames(results_nctree) <- c("correct_main",
                              "correct_composite",
                              "correct_spil")
write.csv(results_nctree, file = "two_main_spillover_four_effects_singular_spil_3000.csv")
