#---------------------------------------------------
# Other.R: Other useful functions
#---------------------------------------------------



#------------------------------------------------------------------
# genmultnet: Generates a simulated clustered network environment

#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Generates a simulated clustered network environment


#' @param  m Number of clusters
#' @param  N Number of units 
#' @param  method method to generate the networks: "ergm" (Exponential Random Graph Models) , "er" (Erdos Renyi) , 
#'       "sf" (Barabasi-Albert model)
#' @param  param If method "er", probability of the ER model
#' @param  varhom Variable to account for homophily 
#' @param  coefergm If method "ergm", coefficients of the ERGM model

#' @return: An adjacency matrix which describes a clustered network environment


genmultnet=function(m,N,method,param,varhom,coefergm){
  comba<-matrix(0,N,N)
  gsize<-N/m
  
  for (i in 1:m){
    
    if(method=="ergm"){
      test.net <- network(gsize, directed = FALSE, density = 0)
      test.net%v%"x1" = varhom[(gsize*i-(gsize-1))]
      g.sim <- simulate(test.net ~ nodematch("x1") + edges,
                        coef = coefergm)
      comba[(gsize*i-(gsize-1)):(gsize*i),(gsize*i-(gsize-1)):(gsize*i)]<-as.matrix(g.sim) 
    }
    
    if(method=="er" | method=="sf"){
      if(method=="er"){
        g=erdos.renyi.game(gsize, p=param, type = "gnp")}
      if(method=="sf"){
        g=barabasi.game(gsize)}
      adj<-as.matrix(get.adjacency(g))
      comba[(gsize*i-(gsize-1)):(gsize*i),(gsize*i-(gsize-1)):(gsize*i)]<-adj
    }
    
  }
  
  return(comba)
}

#---------------------------------------------------
# expand.grid.unique: generates a grid, with no equal elements


#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Create a data frame from all combinations of the supplied vectors or factors,
#' while omitting rows characterized by equal elements

#' @param  x A vector
#' @param  y A vector
#' @param  include.equals=FALSE


#' @return A data frame with all combinations of the elements of the vectors, 
#' with no rows characterized by equal elements


expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#---------------------------------------------------
# sharedn: computes the number of shared neighbors


#' @title
#' Heterogeneous Treatment and Spillover Effects 
#' under Clustered Network Interference

#' @description 
#' Computes the number of shared neighbors between unit i and unit j

#' @param  i Unit ID
#' @param  j Unit ID
#' @param  Nel List of N elements - where N is the sample size - 
#' where each element i contains the IDs of the direct neighbors of unit i


#' @return A numeric value representing the number of shared neighbors between unit i and j


sharedn=function(i,j,Nel){
  return(length(intersect(Nel[[i]],Nel[[j]])))
}
