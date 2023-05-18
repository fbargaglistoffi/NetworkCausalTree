#----------------------------------------------------------
# DGP.R: Function to generate simulated data to run the NCT
#------------------------------------------------------------------

#---------------------------------------------------
# NCT_data_generating_process: returns a Network Causal Tree

#' @title
#' Data Generating Process for Network Causal Tree
#'
#'
#' @description
#' Generates data to be used to generate a Network Causal Tree
#'



#' @param  N Sample size
#' @param  m Number of clusters
#' @param  K Number of binary regressors
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param taui Size of the tretment effects 1000 and 1101
#' @param het If het=TRUE the effects 1000 and 1101 are heterogeneous with respect to the first regressor
#' so that the effects equal +taui in the subpopulation such that X1=0 and -taui in the subpopulation such that X1=1. If het=FALSE you do not introduce heterogeneity and
#' the effects  equal +taui in the whole population.
#' @param  method_networks method to generate the m networks: "ergm" (Exponential Random Graph Models) , "er" (Erdos Renyi) ,
#'       "sf" (Barabasi-Albert model)
#' Note: in this function, clusters have the same size, so N should be a multiple of m
#' @param  param_er If method "er", probability of the ER model
#' @param  var_homophily_ergm Variable to account for homophily in the ERGM model
#' @param  coef_ergm If method "ergm", coefficients of the ERGM model


#' @return A list composed by three elements: i) the first element is the N x N Adjacency matrix;
#' ii) the second element is the dataset N x 6  where the columns provide information
#' on the individual intervention W, on the neighborhood intervention G, on the outcome Y,
#' on the group membership M, #' on the degree Ne, on the probability to be assigned to the active individual intervention;
#' iii) the third element is the N x K covariates matrix X

NCT_data_generating_process=function(N,m,K=5,p,het=TRUE,taui=2,method_networks="er",
                                     param_er=0.1,
                                     coef_ergm=NULL,
                                     var_homophily_ergm=NULL){



#generate covariates
X=NULL
for(k in 1:K){
x=rbinom(N,1,0.5)
X=cbind(X,x)
colnames(X)[k]<-paste0(colnames(X)[k],k)
}

#Generate m networks
gsize=N/m
simnet<-genmultnet(N=N,m=m,method_networks=method_networks,
                   param_er = param_er,
                   coef_ergm = coef_ergm,
                   var_homophily_ergm = var_homophily_ergm)

#group Indicator
M=c(rep(1:m,gsize))
M=sort(M)
levels(M)<-c(1:m)
net<-graph_from_adjacency_matrix(simnet)



#randomly assign unit to treatment arms
treat=c()
for( i in 1:N){
  treat[i]<-rbinom(1,1,prob=p[i])
}


#take the whole adiac_matrix
adiac_matrix<-simnet

neigh<-rowSums(adiac_matrix)
#Compute number of treated neigh and consequently Gi
num_tr_neigh <-as.vector(adiac_matrix %*% treat)
neightreat=rep(1,N) #Gi
neightreat[which(num_tr_neigh==0)]<-0

#Pass to the standard notation
w<-treat[neigh>0]
g<-neightreat[neigh>0]
M<-M[neigh>0]
X<-X[neigh>0,]
p<-p[neigh>0]
N<-length(w)
adiac_matrix<-adiac_matrix[neigh>0,neigh>0]
neigh_red<-neigh[neigh>0]



if(het==TRUE){

x1<-X[,1]
tau = rep(0, N)
tau[x1==0] = taui
tau[x1==1] = - taui

## Generate Treatment Effects
y0 = rnorm(N,sd=0.01)
y1 = y0 + tau
## Generate Outcome
y = y0*(1-w) + y1*w

}

if(het==FALSE){

  tau=rep(taui, N)
  ## Generate Treatment Effects
  y0 = rnorm(N,sd=0.01)
  y1 = y0 + tau
  ## Generate Outcome
  y = y0*(1-w) +  y1*w
}



NCT_dgp<-list()
NCT_dgp[[1]]<-adiac_matrix
NCT_dgp[[2]]<-data.frame(
  W=w,
  G=g,
  Y=y,
  M=M,
  Ne=neigh_red,
  p=p
)
NCT_dgp[[3]]<-X

return(NCT_dgp)

}
