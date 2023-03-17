#----------------------------------------------------------
# Variances.R: Functions to compute the estimated variances of the effects
# of interest
#------------------------------------------------------------------
#------------------------------------------------------------------------
# Vary: computes the estimated variance of the average potential outcomes

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the estimated variance of the average potential outcome 
#' related to a given level of the joint intervention
#' 
#' @param  N Sample size 
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  W N x 1 vector, Individual Treatment 
#' @param  G N x 1 vector, Neighborhood Treatment 
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree  
#' @param  Nel List of N elements - where N is the sample size - 
#' where each element i contains the IDs of the direct neighbors of unit i



#' @return A numeric value corresponding to the estimated variance of the 
#'  average potential outcome under the level w,g of the joint intervention





Vary=function(N,w,g,Y,W,G,p,Ne,Nel){
  varzero=NULL
  variab=c()
  if (length(which(W==w & G==g))>1){
    pairs<-expand.grid.unique(which(W==w & G==g),which(W==w & G==g),include.equals = FALSE)
    for (k in 1:nrow(pairs)){
      i=pairs[k,1]
      j=pairs[k,2]
      if ( sharedn(i=i,j=j,Nel=Nel)>0 ){
        variab=c(varzero,sum(((pij(i,j,w,g,w,g,Ne,Nel,p=p)-pi(i,w,g,p,Ne)*pi(j,w,g,p,Ne))/
                                (pij(i,j,w,g,w,g,Ne,Nel,p=p)))*
                               (Y[i]/pi(i,w,g,p,Ne))*(Y[j]/pi(j,w,g,p,Ne))))
      }}
    vary=sum((1-pi(which(W==w & G==g),w,g,p,Ne))*
               (Y[which(W==w & G==g)]/pi(which(W==w & G==g),w,g,p,Ne))^2) +  
      sum(variab)
  }else{vary=NA}
  
  return(vary)
}  


#------------------------------------------------------------------------
# Covy: computes the estimated covariance of the average potential outcomes

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the estimated covariance of the average potential outcome 
#' related to two given levels of the joint intervention
#' 
#' @param  N Sample size 
#' @param  w1 Individual Treatment level - 1
#' @param  g1 Neighborhood Treatment level - 1
#' @param  w2 Individual Treatment level - 2
#' @param  g2 Neighborhood Treatment level -2
#' @param  W N x 1 vector, Individual Treatment 
#' @param  G N x 1 vector, Neighborhood Treatment 
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree  
#' @param  Nel List of N elements - where N is the sample size - 
#' where each element i contains the IDs of the direct neighbors of unit i



#' @return A numeric value corresponding to the estimated covariance of the 
#'  average potential outcome under the level w1,g1 and w2,g2 of the joint intervention



Covy=function(w1,g1,w2,g2,N,Y,W,G,p,Ne,Nel){
  varzero=NULL
  variab=c()    
  if(length(which(W==w1 & G==g1))>1 & length(which(W==w2 & G==g2))>1){
    pairs<-expand.grid.unique(which(W==w1 & G==g1),which(W==w2 & G==g2),include.equals = FALSE)
    for (k in 1:nrow(pairs)) {
      i=pairs[k,1]
      j=pairs[k,2]
      if(sharedn(i,j,Nel=Nel)>0){
        variab=c(varzero,sum(1/pij(i,j,w1,g1,w2,g2,Ne,Nel,p=p)*
                               (Y[i]/pi(i,w1,g1,p,Ne))*(Y[j]/pi(j,w2,g2,p,Ne)))*
                   (pij(i,j,w1,g1,w2,g2,Ne,Nel,p=p)-pi(i,w1,g1,p,Ne)*pi(j,w2,g2,p,Ne)))
        
      }}
    covy=sum(variab)-(sum(((Y[which(W==w1 & G==g1)])^2)/(2*pi(which(W==w1 & G==g1),w1,g1,p,Ne))) +
                        sum(((Y[which(W==w2 & G==g2)])^2)/(2*pi(which(W==w2 & G==g2),w2,g2,p,Ne))))  
  }else{covy=NA}
}


#---------------------------------------------------
# Vartau1000, Vartau1101 ,Vartau0100 ,Vartau1110 : 
# Functions to compute the estimated variance of the
# four effects of interest 

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the estimated variance of the effect of interest 
#' 
#' @param  N Sample size 
#' @param  W N x 1 vector, Individual Treatment 
#' @param  G N x 1 vector, Neighborhood Treatment 
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree  
#' @param  Nel List of N elements - where N is the sample size - 
#' where each element i contains the IDs of the direct neighbors of unit i



#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.




Vartau1000=function(N,Y,W,G,p,Ne,Nel){
  
  vary10<-Vary(w=1,g=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  vary00<-Vary(w=0,g=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  covy10y00 <-Covy(w1=1,g1=0,w2=0,g2=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  if(any(is.na(c(vary10,vary00)))){
    var1000=NA
  }else{ 
    var1000=abs((1/(N^2))*(vary10+vary00-2*covy10y00))
  }
  return(var1000) 
}

Vartau1101=function(N,Y,W,G,p,Ne,Nel){
  
  vary11<-Vary(w=1,g=1,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  vary01<-Vary(w=0,g=1,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  covy11y01 <-Covy(w1=1,g1=1,w2=0,g2=1,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  if(any(is.na(c(vary11,vary01)))){
    var1101=NA
  }else{  
    var1101=abs((1/(N^2))*(vary11+vary01-2*covy11y01))
  }
  return(var1101) 
}

Vartau1110=function(N,Y,W,G,p,Ne,Nel){
  
  vary11<-Vary(w=1,g=1,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  vary10<-Vary(w=1,g=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  covy11y10 <-Covy(w1=1,g1=1,w2=1,g2=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  if(any(is.na(c(vary11,vary10)))){
    var1110=NA
  }else{ 
    var1110=abs((1/(N^2))*(vary11+vary10-2*covy11y10))
  }
  return(var1110) 
}

Vartau0100=function(N,Y,W,G,p,Ne,Nel){
  
  vary01<-Vary(w=0,g=1,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  vary00<-Vary(w=0,g=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  covy01y00 <-Covy(w1=0,g1=1,w2=0,g2=0,N=N,W=W,G=G,p=p,Ne=Ne,Nel=Nel,Y=Y)
  
  if(any(is.na(c(vary01,vary00)))){
    var0100=NA
  }else{
    var0100=abs((1/(N^(2)))*(vary01+vary00-2*covy01y00))
  }
  return(var0100) 
}