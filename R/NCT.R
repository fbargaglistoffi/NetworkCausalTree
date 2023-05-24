#' @title
#' Network Causal Tree

#' @description
#' Returns a Network Causal Tree, with the corresponding estimates
#'
#' @param X N x K Observed Covariates Matrix.
#' @param Y N x 1 Observed Outcome vector.
#' @param W N x 1 Individual Treatment vector.
#' @param effweights Treatment Effect weights vector (4 elements):
#' - alpha: weight associated to the treatment effect 1000 (effect of the
#' individual treatment, with the neighborhood treatment set at 0),
#' - beta: weight associated to the treatment effect 1101 (effect of the
#' individual treatment, with the neighborhood treatment set at 1),
#' - gamma: weight associated to the spillover effect 1110 (effect of the
#' neighborhood treatment, with the individual treatment set at 1),
#' - delta: weight associated to the spillover effect 0100 (effect of the
#' neighborhood treatment, with the individual treatment set at 0).
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
#' @return A data.frame describing the obtained Network Causal Trees.
#' Each row represents a partition (of a specific tree) with 10/14 entries.
#' Columns summary:
#' - `OF`: value of the OF in the corresponding partition,
#' - `NOBS_TR`: number of training observations in the partition,
#' - `FILTER`: values of the covariates `X` that identify the partition,
#' - `IDTREE`: tree ID,
#' - `NUMTREE`: number of trees identifying the partition,
#' - `NOBS_EST`: number of estimation observations in the partition,
#' - `EFF1000_EST`: estimated 1000 effects in the partitions,
#' - `EFF1101_EST`: estimated 1101 effects in the partitions,
#' - `EFF1110_EST`: estimated 1110 effects in the partitions,
#' - `EFF0100_EST`: estimated 0100 effects in the partitions.
#' Additional columns summary (only if output = "Estimation"):
#' - `SETAU1000_EST`: estimated std. error of the 1000 effect in the partition,
#' - `SETAU1101_EST`: estimated std. error of the 1101 effect in the partition,
#' - `SETAU1110_EST`: estimated std. error of the 1110 effect in the partition,
#' - `SETAU0100_EST`: estimated std. error of the 0100 effect in the partition.
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
