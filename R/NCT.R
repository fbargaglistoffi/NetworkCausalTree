#' @title
#' Network Causal Tree

#' @description
#' Returns a Network Causal Tree, with the corresponding estimates
#'
#' @param X N x K Observed Covariates Matrix.
#' @param Y N x 1 Observed Outcome vector.
#' @param W N x 1 Individual Treatment vector.
#' @param effweights Vector including the 4 effect weight:
#' - alpha weight associated to the effect 1000,
#' - beta weight associated to the effect 1101,
#' - gamma weight associated to the effect 1110,
#' - delta weight associated to the effect 0100.
#' @param A N x N Adjacency matrix.
#' @param G N x 1 Neighborhood Treatment vector.
#' @param M N x 1 Cluster Membership vector.
#' @param p  N x 1 Probability to be assigned to the active individual
#' intervention vector.
#' @param mdisc Number of clusters to be assigned to the discovery set.
#' @param mest Number of clusters to be assigned to the estimation set only.
#' @param minpopfrac Ratio of the discovery set population to be included while
#' sprouting the tree.
#' @param fracpredictors Quote of the predictors to be included while sprouting
#' the tree
#' @param n_trees Number of Trees.
#' @param minsize Minimum number of observaztions for each level of the joint
#' intervention
#' to be required in the leafs
#' @param depth Depth of the tree.
#' @param method Method to compute the objective function: "singular" for NCT
#' targeted to one single effect; "composite" for NCT targeted to multiple
#' effects; "penalized" for a OF computed while considering a single effect only
#' and including a penalization term related to the variance.
#' @param output Desired output of the analysis. if output = "detection" only
#' point estimates are computed, if output = "estimation" both estimated effects
#' and variances are computed
#'
#' @return A Network Causal Tree - NCT - object. An NCT object is a data frame reporting the
#' results of the Network Causal Trees, where each tree is characterized by a set of partitions.
#' Specifically the GOF variable includes the value of the GOF in the corresponding partition, the NOBS_TR column
#' includes the number of obs belonging to the training set and to the given partition,
#' the column FILTER includes the values of the Xs that identify the given partition, the
#' IDTREE variable identifies the id of the tree ,  the
#' NUMTREE variable represents the number of trees that identify the given partition,
#' the NOBS_EST column includes the number of obs belonging
#' to the estimation set and to the given partition, the columns EFF1000_EST, EFF1101_EST,  EFF1110_EST,
#' EFF0100_EST reporting the estimated effects in the partitions. If output = "Estimation" the data frame
#' also include the columns SETAU1000_EST, SETAU1101_EST,  SETAU1110_EST, SETAU0100_EST reporting the estimated standard errors of the effects in the partitions
#'
#' @import stringi
#' @import igraph
#' @import statnet
#' @import ergm
#' @import plyr
#' @import data.tree
#' @import stats
#'
#' @export
#'
NetworkCausalTrees=function(X, Y, W,
                            effweights = c(0.25,0.25,0.25,0.25),
                            A = NULL,
                            G = NULL,
                            M = NULL,
                            p = NULL,
                            mdisc = 25,
                            mest = 15,
                            minpopfrac = 1,
                            fracpredictors = 1,
                            n_trees = 1,
                            depth = 3,
                            minsize = 10,
                            method = "singular",
                            output = "estimation"){

  N <- length(W)
  m <- length(unique(M))

  Ne <- rowSums(A)

  # get input weights
  alpha <- effweights[1]
  beta <- effweights[2]
  gamma <-effweights[3]
  delta <- effweights[4]


  #error messages
  if(alpha+beta+gamma+delta!=1){
    stop('weights have to sum up to one')
  }

  if((length(which(effweights>0))>1) & (method=="singular" | method=="penalized")){
    stop('if method is set to singular only one effect should have positive weight')
  }

  if(1 %in% effweights & method=="composite"){
    stop('composite gof is computed if at least two effects are investigated')
  }

  if(is.null(G) & is.null(A)){
    stop('independently of the output, you have to specify either the G vector or the Adiacency Matrix')
  }

  if(is.null(Ne) & is.null(A)){
    stop('independently of the output, you have to specify either the Degree vector or the Adiacency Matrix')
  }

  if(output=="estimation" & is.null(A)){
    stop('if output is set estimation you have to specify the adiacency matrix')
  }

  if(is.null(Ne)){
    Ne<-rowSums(A)
  }


  if(is.null(G)){
    nt <-as.vector(A %*% W)
    G=rep(1,N)
    G[nt==0]<-0}

  Nel<-vector(mode = "list", length = N)
  for( i in  1:N){
    Nel[[i]]<-which(A[i,]>0)
  }

  if(is.null(Ne)){
    stop('if you dont input the Adiacency Matrix you must specify the degree vector')
  }



  data <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X,M=M)

  Peff<-popeff(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)

  sampgroup_train=sample(1:m,size=mdisc,replace=FALSE)

  trees <- plyr::raply(
    n_trees,
    sproutnetctree(method=method,sampgroup=sampgroup_train,fracpredictors=fracpredictors,m=m,minpopfrac=minpopfrac,
                   depth=depth,minsize=minsize,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                   N=N,W=W,G=G,Y=Y,X=X,M=M,Ne=Ne,p=p,Peff=Peff,Nel=Nel),
    .progress = "text" )

  Forest<-data.frame(IDTREE=NA,NODE = NA, GOF=NA, NOBS = NA, FILTER = NA, TERMINAL = NA,
                     stringsAsFactors = FALSE)


  for(i in 1:n_trees){
    tree=NA
    tree<-cbind(rep(i,nrow(as.data.frame(trees[i]))),as.data.frame(trees[i]))
    colnames(tree)<-colnames(Forest)
    Forest<-rbind(Forest,tree)
  }

  Forest=Forest[-1,]

  Results<-data.frame(GOF=Forest$GOF,NOBS=Forest$NOBS ,  FILTER =  c("NA",as.vector(na.omit(unique(Forest$FILTER)))),
                      NUMTREE=NA, IDTREE=Forest$IDTREE, stringsAsFactors = FALSE)

  for(j in as.vector(na.omit(unique(Forest$FILTER)))){
    Results$NUMTREE[Results$FILTER==j]<-length(unique(Forest$IDTREE[which(Forest$FILTER==j)]))
    Results$NUMTREE[Results$FILTER=="NA"]=length(unique(Forest$IDTREE[is.na(Forest$FILTER)]))
  }

  if(nrow(Results)==1){
    warning('No split has been made')
  }

  Results_est<-alleffect(output=output,tree_info = Results,N=N,W=W,G=G,Y=Y,X=X,Ne=Ne,p=p,Nel=Nel,minsize=minsize)
  return(Results_est)
}
