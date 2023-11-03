#' @title
#' Estimated Variance of the Average Potential Outcomes
#'
#' @description
#' Computes the estimated variance of the average potential outcome
#' related to a given level of the joint intervention - denoted with (w,g) -.
#'
#' @param  N Sample size
#' @param  w Individual Treatment level
#' @param  g Neighborhood Treatment level
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the
#'  average potential outcome under the level w,g of the joint intervention


Vary = function(N, w, g, Y, W, G, 
                p, Ne, Ne_list) {
  
  varzero = NULL
  second_element = c()
  
  if (length(which(W == w & G == g)) > 1) {
    
    pairs <- expand.grid.unique(which(W == w & G == g), 
                                which(W == w & G == g),
                                include.equals = FALSE)
    
    for (k in 1 : nrow(pairs)) { 
      
      i = pairs[k, 1]
      j = pairs[k, 2]
      
      if (shared_neigh(i = i, j = j, Ne_list = Ne_list) > 0) {
        
        second_element = c(varzero, sum(((pij(i, j, w, g, w, g, Ne, Ne_list, p = p)
                         -pi(i, w, g, p, Ne) * pi(j, w, g, p, Ne)) /
                         (pij(i, j, w, g, w, g, Ne, Ne_list, p = p))) *
                         (Y[i] / pi(i, w, g, p, Ne)) * (Y[j] / pi(j, w, g, p,Ne))))
      }
      
      }
    
    vary = sum((1 - pi(which(W == w & G==g), w, g, p, Ne)) *
           (Y[which(W == w & G == g)] / pi(which(W == w & G == g), w, g, p, Ne)) ^ 2) +
           sum(second_element)
    
  } else {vary = NA}

  return(vary)
}

# ------------------------------------------------------------------------------

#' @title
#' Estimated Covariance of Average Potential Outcomes
#'
#' @description
#' Computes the estimated covariance of the average potential outcomes
#' related to two given levels of the joint intervention
#'  - denoted with (w1,g1) and (w2,g2), respectively -.
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
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated covariance of the
#'  average potential outcome under the level w1,g1 and w2,g2 of the joint intervention

Covy = function(w1, g1, w2, g2, 
                N, Y, W, G, p, 
                Ne, Ne_list) {
  
  varzero = NULL
  variab = c()
  
  if (length(which(W == w1 & G == g1)) > 1 & 
      length(which(W == w2 & G == g2)) > 1) {
    
     pairs <- expand.grid.unique(which(W == w1 & G == g1),
                                 which(W == w2 & G == g2),
                                 include.equals = FALSE)
     
    for (k in 1:nrow(pairs)) {
      
      i=pairs[k,1]
      j=pairs[k,2]
      
      if (shared_neigh(i, j, Ne_list = Ne_list) > 0) {
        
        variab = c(varzero,
                   sum(1 / pij(i, j, w1, g1, w2, g2, Ne, Ne_list, p = p)*
                      (Y[i] / pi(i, w1, g1, p, Ne)) * 
                      (Y[j] / pi(j, w2, g2, p, Ne))) *
                      (pij(i, j, w1, g1, w2, g2, Ne, Ne_list, p = p) -
                      pi(i, w1, g1, p, Ne) * pi(j, w2, g2, p,Ne)))

      }
      
    }
    
     second_element <- sum(((Y[which(W == w1 & G == g1)]) ^ 2) / (2 * 
                       pi(which(W == w1 & G == g1), w1, g1, p, Ne))) +
                       sum(((Y[which(W == w2 & G == g2)]) ^ 2) / (2 * 
                       pi(which(W == w2 & G == g2), w2, g2, p, Ne))) 
   
     
    covy = sum(variab) - second_element
      
  } else {covy=NA}
}


#-------------------------------------------------------------------------------

#' @title
#' Estimated Variance for the effect 1000
#'
#' @description
#' Computes the estimated variance of the effect 1000. This effect
#'  compares individuals with the joint treatment exposure (w=1,g=0)
#' with those with (w=0,g=0)
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.
#'
Vartau1000 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {

# Variance of the Average Potential Outcome Y(w,g)=Y(1,0)
  vary10 <- Vary(w = 1, g = 0, N = N, W = W, G = G, p = p,
                 Ne = Ne, Ne_list = Ne_list,Y = Y)

# Variance of the Average Potential Outcome Y(w,g)=Y(0,0)    
  vary00 <- Vary(w = 0, g = 0, N = N, W = W, G = G, p = p,
                 Ne = Ne, Ne_list = Ne_list,Y = Y)

# Covariance of the Average Potential Outcomes Y(w,g)=Y(1,0) and Y(w,g)=Y(0,0) 
  covy10y00 <- Covy(w1 = 1, g1 = 0, w2 = 0, g2 = 0,
                    N = N, W = W, G = G, p = p,
                    Ne = Ne, Ne_list = Ne_list, Y = Y)

  if (any(is.na(c(vary10, vary00)))) {
    var1000 = NA
  } else {
    var1000 = abs((1 / (N ^ 2)) * 
                 (vary10 + vary00 - 2 * covy10y00))
  }
  return(var1000)
}


#-------------------------------------------------------------------------------

#' @title
#' Estimated Variance for the effect 1101
#'
#' @description
#' Computes the estimated variance of the effect 1000. This effect
#'  compares individuals with the joint treatment exposure (w=1,g=1)
#' with those with (w=0,g=1)
#'
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.
#'
Vartau1101 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {

  # Variance of the Average Potential Outcome Y(w,g)=Y(1,1)  
  vary11 <- Vary(w = 1, g = 1, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(0,1)
  vary01 <- Vary(w = 0, g = 1, N = N, W = W, G = G,
                 p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)

  # Covariance of the Average Potential Outcomes Y(w,g)=Y(1,1) and Y(w,g)=Y(0,1) 
    covy11y01 <- Covy(w1 = 1, g1 = 1,w2 = 0,g2 = 1,
                      N = N,W = W,G = G,
                      p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)

  if (any(is.na(c(vary11, vary01)))) {
    var1101 = NA
  } else {
    var1101 = abs((1 / (N ^2)) * 
                  (vary11 + vary01 - 2 * covy11y01))
  }
  return(var1101)
}


#-------------------------------------------------------------------------------

#' @title
#' Estimated Variance for the effect 1110
#'
#' @description
#' Computes the estimated variance of the effect 1110. This effect
#'  compares individuals with the joint treatment exposure (w=1,g=1)
#' with those with (w=1,g=0)
#'
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.
#'
Vartau1110 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {
 
  # Variance of the Average Potential Outcome Y(w,g)=Y(1,1)   
  vary11 <- Vary(w = 1, g = 1, N = N, W = W, G = G,
               p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(1,0)
  vary10 <- Vary(w = 1, g = 0, N = N, W = W, G = G,
               p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)

  # Covariance of the Average Potential Outcomes Y(w,g)=Y(1,1) and Y(w,g)=Y(1,0) 
  covy11y10 <- Covy(w1 = 1, g1 = 1, w2 = 1, g2 = 0,
                    N = N, W = W, G = G, p = p,
                    Ne = Ne, Ne_list = Ne_list, Y = Y)

  if (any(is.na(c(vary11, vary10)))) { 
    var1110 = NA
  } else { 
    var1110 = abs((1 / (N ^ 2)) * 
                    (vary11 + vary10 - 2 * covy11y10))
  }
  return(var1110)
}


#-------------------------------------------------------------------------------

#' @title
#' Estimated Variance for the effect 0100
#'
#' @description
#' Computes the estimated variance of the effect 0100. This effect
#'  compares individuals with the joint treatment exposure (w=0,g=1)
#' with those with (w=0,g=0)
#'
#'
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  p N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Ne_list List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#'
#' @return A numeric value corresponding to the estimated variance of the effect
#' of interest.
#'
Vartau0100 = function(N, Y, W, G, 
                      p, Ne, Ne_list) {
  
  # Variance of the Average Potential Outcome Y(w,g)=Y(0,1)
  vary01 <- Vary(w = 0, g = 1, N = N, W = W, G = G,
               p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)

  # Variance of the Average Potential Outcome Y(w,g)=Y(0,0)
  vary00 <- Vary(w = 0, g = 0, N = N, W = W, G = G,
               p = p, Ne = Ne, Ne_list = Ne_list, Y = Y)

  # Covariance of the Average Potential Outcomes Y(w,g)=Y(0,1) and Y(w,g)=Y(0,0) 
   covy01y00 <- Covy(w1 = 0, g1 = 1, w2 = 0, g2 = 0,
                    N = N, W = W, G = G, p = p,
                    Ne = Ne, Ne_list = Ne_list, Y = Y)

  if (any(is.na(c(vary01, vary00)))) {
    var0100 = NA
  } else {
    var0100 = abs((1 / (N ^ (2))) *
                    (vary01 + vary00 - 2 * covy01y00))
  }
  return(var0100)
}
