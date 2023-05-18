#' @title
#' Estimated Effect 1000
#'
#' @description
#' Computes the estimates of the effect of interest
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau1000=function(N,W,G,Y,p,Ne){
  tau1000=1/N*(sum(Y[W==1 & G==0]/pi(which(W==1 & G==0),1,0,p,Ne))-
                 sum(Y[W==0 & G==0]/pi(which(W==0 & G==0),0,0,p,Ne)))
  return(tau1000)
}

#-------------------------------------------------------------------------------

#' @title
#' Estimated Effect 1101

#' @description
#' Computes the estimates of the effect of interest
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau1101=function(N,W,G,Y,p,Ne){
  tau1101=1/N*(sum(Y[W==1 & G==1]/pi(which(W==1 & G==1),1,1,p,Ne))-
              sum(Y[W==0 & G==1]/pi(which(W==0 & G==1),0,1,p,Ne)))
  return(tau1101)
}

#-------------------------------------------------------------------------------

#' @title
#' Estimated Effect 1110
#'
#' @description
#' Computes the estimates of the effect of interest
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau1110=function(N,W,G,Y,p,Ne){
  tau1110=1/N*(sum(Y[W==1 & G==1]/pi(which(W==1 & G==1),1,1,p,Ne))-
                 sum(Y[W==1 & G==0]/pi(which(W==1 & G==0),1,0,p,Ne)))
  return(tau1110)
}

#-------------------------------------------------------------------------------

#' @title
#' Estimated Effect 0100
#'
#' @description
#' Computes the estimates of the effect of interest
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric value corresponding to the estimate of the effect
#' of interest.
#'
EffTau0100=function(N,W,G,Y,p,Ne){
  tau0100=1/N*(sum(Y[W==0 & G==1]/pi(which(W==0 & G==1),0,1,p,Ne))-
                 sum(Y[W==0 & G==0]/pi(which(W==0 & G==0),0,0,p,Ne)))
  return(tau0100)
}

#-------------------------------------------------------------------------------

#' @title
#' Estimated Effects in the whole sample
#'
#' @description
#' Computes the four Estimated Effects in the Whole Population
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#'
#' @return A numeric vector made up by four elements, representing the
#' estimated effects in the whole population: the first element is the effect 1000,
#' the second one refers to the effect 1101, the third one refers to 1110 and the
#' fourth refers to 0100
#'
popeff=function(N,W,G,Y,p,Ne){
  PEffTau1000<-EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  PEffTau1101<-EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  PEffTau1110<-EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  PEffTau0100<-EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  Peff<-c(PEffTau1000,PEffTau1101,PEffTau1110,PEffTau0100)
  return(Peff)
}
