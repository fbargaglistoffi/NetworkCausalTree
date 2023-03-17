#----------------------------------------------------------
# Probabilities.R: Functions to compute marginal and joint probabilities
#------------------------------------------------------------------



#---------------------------------------------------
# pi: Individual marginal probabilities 

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the individual marginal probability to be exposed to a given level of the joint intervention

#' @param  i Unit ID
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree  



#' @return A numeric value bounded between 0 and 1 measuring 
#' the individual marginal probability to be exposed 
#' to the level w,g of  the joint intervention



pi<-function(i,w,g,p,Ne){
  pi<-(p[i]^w)*
    (1-p[i])^(1-w)*
    ((1-(1-p[i])^Ne[i])^g)*
    ((1-p[i])^Ne[i])^(1-g)
  return(pi)
}

#---------------------------------------------------
# pij: Individual joint probabilities 

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the individual joint probability of a given pair of units
#'  to be simultaneously exposed to certain levels of the joint intervention

#' @param  i Unit ID
#' @param  j Unit ID
#' @param  wi Individual Treatment level
#' @param  wj Individual Treatment level
#' @param  gi Neighborhood Treatment level
#' @param  gj Neighborhood Treatment level
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree   
#' @param  Nel List of N elements - where N is the sample size - 
#' where each element i contains the IDs of the direct neighbors of unit i




#' @return A numeric value bounded between 0 and 1 measuring 
#' the individual joint probability of units i and j to be simultaneously
#' exposed to the level wi,gi and wj,gj of  the joint intervention, respectively

pij<-function(i,j,wi,wj,gi,gj,Ne,Nel,p){
  pij=  ((p[i]^wi)*(1-p[i])^(1-wi))*
    ((p[j]^wj)*(1-p[j])^(1-wj))*
    ((1-(1-p[i])^(Ne[i]-ifelse(Ne[i]>=Ne[j] & Ne[i]>sharedn(i,j,Nel=Nel) ,sharedn(i,j,Nel=Nel),0)))^gi)*
    ((1-p[i])^(Ne[i]-ifelse(Ne[i]>=Ne[j] & Ne[i]>sharedn(i,j,Nel=Nel),sharedn(i,j,Nel=Nel),0)))^(1-gi)*
    ((1-(1-p[j])^(Ne[j]-ifelse(Ne[i]<Ne[j],sharedn(i,j,Nel=Nel),0)))^gj)*
    ((1-p[j])^(Ne[j]-ifelse(Ne[i]<Ne[j],sharedn(i,j,Nel=Nel),0))^(1-gj))
  
  return(pij)  
  
}
  
  
  
