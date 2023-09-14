#' @title
#' Objective Function (OF)
#'
#' @description
#' Computes the measure of the Objective Function
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param p  N x 1 vector,Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param Peff 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param vartot - to be included if method = "penalized" - whole variance
#' @param leafs - to be included if method = "penalized" - number of leafs
#'
#' @return A numeric value corresponding to the computed  Objective Function
#'
OF_Value=function(method,alpha,beta,gamma,delta,N,W,G,Y,p,Ne,Nel,Peff,vartot,leafs){
  inof=NULL
  if(method=="composite"){
    inof<-alpha*(((EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[1])^2)+
      beta*(((EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[2])^2)+
      gamma*(((EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[3])^2)+
      delta*(((EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[4])^2)
  }

  if(method=="penalized")
  {l=leafs
  inof<-
    alpha*(
      ((EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2) -
        2/l*sum(c(vartot,Vartau1000(N=N,W=W,Y=Y,G=G,p=p,Ne=Ne,Nel=Nel)))
    ) +
    beta*(
      ((EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)-
        2/l*sum(c(vartot,Vartau1101(N=N,W=W,Y=Y,G=G,p=p,Ne=Ne,Nel=Nel)))
    ) *
    gamma*(
      ((EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)-
        2/l*sum(c(vartot,Vartau1110(N=N,W=W,Y=Y,G=G,p=p,Ne=Ne,Nel=Nel)))
    )+
    delta*(
      ((EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)-
        2/l*sum(c(vartot,Vartau0100(N=N,W=W,Y=Y,G=G,p=p,Ne=Ne,Nel=Nel)))
    )
  }

  if (method=="singular")
  {inof<-alpha*((EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
    beta*((EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
    gamma*((EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
    delta*((EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)
  }
  return(inof)
}



#-------------------------------------------------------------------------------

#' @title
#' Split Objective Function
#'
#' @description
#' Splits the sample where the  Objective Function is maximized
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param N Sample size
#' @param W N x 1 vector, Individual Treatment
#' @param G N x 1 vector, Neighborhood Treatment
#' @param Y N x 1 vector, Observed Outcome
#' @param X N x K matrix, Observed Covariate Matrix
#' @param p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param Ne N x 1 vector, Degree
#' @param Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param Peff 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param vartot - to be included if method = "penalized" - whole variance
#' @param leafs - to be included if method = "penalized" - number of leafs
#'
#' @return A numeric vector made up by three elements: the first one identifies the
#' value of the  Objective Function is max , the second one reports the value of the variable that maximizes the  Objective Function,
#' the third one reports the corresponding variable
#'
OF_Split=function(method,alpha,beta,gamma,delta,N,W,G,Y,X,p,Ne,Nel,Peff,vartot,leafs){

  #initialize
  of<- c()
  name <-c()
  values<-c()
  #check all the predictors
  for(j in 1:dim(X)[2]){
    x=X[,j]
    splits <- sort(unique(x))
    valuesx<-c()
    ofx <- c()
    namesx<-c()
    for(i in seq_along(splits)){
      sp <- splits[i]
      if (all(as.numeric(table(W[x>=sp],G[x>=sp]))>2) & all(as.numeric(table(W[x<sp],G[x<sp]))>2)) {
        ofx[i]<- 1/2*(OF_Value(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta,N=length(which(x<sp)),W=W[x < sp],G=G[x < sp],Y=Y[x < sp],Ne=Ne[x < sp],p=p[x < sp],Peff=Peff,Nel=Nel[x < sp],vartot=vartot,leafs=leafs) +
                        OF_Value(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta,N=length(which(x>=sp)),W=W[x >= sp],G=G[x >= sp],Y=Y[x >= sp],Ne=Ne[x >= sp],p=p[x >= sp], Peff = Peff,Nel=Nel[x >= sp],vartot=vartot,leafs=leafs))} else {ofx[i]<-0}
    }
    namex=rep(colnames(X)[j],length(unique(x)))
    valuesx=c(sort(unique(x)))

    #append all the computed values of the OF, all the values that have defined the split and the name of the variable that has been used.
    of=c(of,ofx)
    name=c(name,namex)
    values=c(values,valuesx)
  }

  if(all(is.na(of))){
    ofres<-NA
    splitres<-NA
    varres<-NA
  }else{
    ofres<-max(na.omit(of))
    splitres<-values[which.max(of)]
    varres<-name[which.max(of)]
  }

    return(c(of = ofres , split = splitres , var=varres))

}


#-------------------------------------------------------------------------------

#' @title
#' Identification of Partitions of the Network Causal Tree
#'
#' @description
#' Identifies the partitions of a Network Causal Tree
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  X N x K matrix, Observed Covariate Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  Peff 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize Minimum number of observations for each level of the joint intervention
#' to be required in the leafs
#' @param depth Depth of the tree
#'
#' @return A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the  Objective Function in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' valyes of the Xs to identify the given partition, the column TERMINAL reports the
#' 'role' of the node - parent or leaf -
#'
netctree <- function(method,alpha,beta,gamma,delta,depth,minsize,N,W,G,Y,X,p,Ne,Nel,Peff)
{

  data_tree <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)
  do_splits <- TRUE

  #CREATE OUTPUT DATASET
  tree_info <- data.frame(NODE = 1, OF=0, NOBS = nrow(data_tree), FILTER = NA, TERMINAL = "SPLIT",
                          stringsAsFactors = FALSE)
  vartot=NULL
  while(do_splits){
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")

    for (j in to_calculate) {

      if (!is.na(tree_info[j, "FILTER"])){
        texts=tree_info[j, "FILTER"]
        this_data <- subset(data_tree, eval(parse(text=texts)))
      } else {
        this_data<- data_tree}

      leafs=nrow(tree_info)
      #SPLIT WRT OF OVER THE SUBSET
      splitting <- OF_Split(method=method,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                             N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,
                             X=this_data[, grepl("X.",names(this_data))],
                             Ne=Ne[this_data$idunit],p=p[this_data$idunit],
                             Nel=Nel[this_data$idunit],
                             Peff = Peff,leafs=leafs,vartot = vartot)

      if (any(is.na(splitting))) {
        split_here <- rep(FALSE, 2)
        print('splits has stopped couse OF is all NA')
      }else{

        #GET THE MAX OF
        maxof <- as.numeric(splitting[1])
        mn <- max(tree_info$NODE)


        #PASTE FILTER RULES
        tmp_filter <- c(paste("data_tree$",splitting[3], ">=","(" ,
                              as.numeric(splitting[2]),")",sep=""),
                        paste("data_tree$",splitting[3], "<", "(",
                              as.numeric(splitting[2]),")",sep=""))
      }
      #CHECK IF THE CURRENT SPLITTING RULE HAS ALREADY BEEN USED
      split_here  <- !sapply(tmp_filter,
                             FUN = function(x,y) any(grepl(x, x = y)),
                             y = tree_info$FILTER)

      #APPEND SPLITTING RULES
      if (!is.na(tree_info[j, "FILTER"])) {
        tmp_filter  <- paste(tree_info[j, "FILTER"],
                             tmp_filter, sep = " & ") }

      #COUNT THE NUMBER OF NODES IN THE CURRENT NODE
      tmp_nobs <- sapply(tmp_filter,
                         FUN = function(i, x) {
                           nrow(subset(x = x, subset = eval(parse(text = i))))
                         },
                         x = this_data)


      #STOP IF INSUFFICIENT MINSIZE
      if (any(as.numeric(table(this_data$W,this_data$G))< minsize)) {
        split_here <- rep(FALSE, 2)
        print('split has stopped for insufficient minsize')
      }


      #STOP IF THE THREE IS OVER THE DEPTH
      depth_tree<-as.numeric(stri_count_regex(tree_info[j, "FILTER"], "X."))
      if (depth_tree>=depth & !is.na(depth_tree)){
        split_here <- rep(FALSE, 2)
        print('split has stopped for reached depth')
      }



      #CREATE CHILDREN DATASET
      children <- data.frame(NODE = c(mn+1, mn+2),
                             OF=c(rep(maxof,2)),
                             NOBS = tmp_nobs,
                             FILTER = tmp_filter,
                             TERMINAL = rep("SPLIT", 2),
                             row.names = NULL)[split_here,]

      if(method=="penalized"){
        varchild=alpha*(Vartau1000(N=nrow(this_data),W=this_data$W,Y=this_data$Y,G=this_data$G,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nel[this_data$idunit]))+
          beta*(Vartau1101(N=nrow(this_data),W=this_data$W,Y=this_data$Y,G=this_data$G,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nel[this_data$idunit]))+
          gamma*(Vartau1110(N=nrow(this_data),W=this_data$W,Y=this_data$Y,G=this_data$G,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nel[this_data$idunit]))+
          delta*(Vartau0100(N=nrow(this_data),W=this_data$W,Y=this_data$Y,G=this_data$G,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nel[this_data$idunit]))
        vartot=c(vartot,varchild)
      }

      #ASSIGN PARENTS OR LEAF IF SPLIT_HERE IS ACTIVE OR NOT
      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")
      tree_info <- rbind(tree_info, children)
      #STOP SPLITTING IF THERE AREN'T NODES WITH "SPLIT"
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    }

  }

  return(tree = tree_info)
}


#-------------------------------------------------------------------------------

#' @title
#' Generation of the Network Causal Tree
#'
#' @description
#' Sprouts the network causal tree, eventually including a fraction of the initial -discovery-
#' sample or a fraction of predictors.
#'
#' @param method method to compute the Objective function: "singular" for NCT targeted to one single effect;
#' "composite" for NCT targeted to multiple effects; "penalized" for a OF computed while
#' considering a single effect only and including a penalization term related to the variance
#' @param alpha weight associated to the effect 1000
#' @param beta weight associated to the effect 1101
#' @param gamma weight associated to the effect 1110
#' @param delta weight associated to the effect 0100
#' @param  N Sample size
#' @param sampgroup Clusters assigned to the discovery set
#' @param m Total number of clusters
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  M N x 1 vector, Cluster Membership
#' @param  X N x K matrix, Observed Covariates Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  Peff 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize minimum number of observaztions for each level of the joint intervention
#' to be required in the leafs
#' @param depth depth of the tree
#'
#' @return A list containing two elements 1- A data frame describing the Network Causal Tree. Specifically,
#' the data frame includes the NODE column identifying the node number, the OF variable
#' including the value of the OF in the corresponding partition, the NOBS column
#' including the number of obs in the given partition, the column FILTER including the
#' valyes of the Xs to identify the given partition,, the column TERMINAL reports the
#' 'role' of the node - parent or leaf -. 2 - a vector that includes the
#' predictors employed while sprouting the tree
#'
sproutnetctree=function(method,minpopfrac,fracpredictors,sampgroup,m,alpha,beta,gamma,delta,depth,minsize,N,W,G,Y,X,M,p,Ne,Peff,Nel){

  # coerce to data.frame
  data <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X,M=M)

  #SAMPLE PREDICTORS
  samppredictors <- sort(sample(which(grepl("X.",names(data))),
                                size = ceiling(length(which(grepl("X.",names(data)))) * 1),
                                replace = FALSE))

  #CREATE SUBSET OF DATA ABOVE WHICH IS USED TO BUILD UP THE TREE AND BUILD THE TREE
  datasample<-data[which(M %in% sampgroup),c(1:4,samppredictors,ncol(data))]
  sampunit<- nrow(datasample)

  datasample<-datasample[sample(1:nrow(datasample), size = sampunit, replace = FALSE),]
  datasample<-datasample[order(datasample$idunit) ,]
  sampleid<-unique(datasample$idunit)
  N=length(sampleid)
  W=as.numeric(datasample$W)
  G=as.numeric(datasample$G)
  Y=as.numeric(datasample$Y)
  X=as.matrix(datasample[,-c(1:4,dim(datasample)[2])])
  colnames(X)=sub("X.","",colnames(X))

  tree_info<-netctree(method=method,alpha=alpha,beta=beta,gamma=gamma,delta=delta,N=length(sampleid),
                      depth=depth,minsize=minsize,W=W,G=G,Y=Y,X=X,
                      p=p[sampleid],Ne=Ne[sampleid],Nel=Nel[sampleid],Peff=Peff)

  return(list(tree=tree_info,predictors_used=colnames(data[,samppredictors])))
}

#-------------------------------------------------------------------------------
#' @title
#' Computation of the Effects in all NCT partitions
#'
#' @description
#' Computes the estimates in all the partitions identified by the Network Causal Tree
#'
#' @param output Desired output of the analysis. if output = "detection" only point estimates
#' are computed, if output = "estimation" both estimated effects and variances are computed
#' @param tree_info An NCT data frame
#' @param  N Sample size
#' @param  W N x 1 vector, Individual Treatment
#' @param  G N x 1 vector, Neighborhood Treatment
#' @param  Y N x 1 vector, Observed Outcome
#' @param  X N x K matrix, Observed Covariates Matrix
#' @param  p  N x 1 vector, Probability to be assigned to the active individual intervention
#' @param  Ne N x 1 vector, Degree
#' @param  Nel List of N elements - where N is the sample size -
#' where each element i contains the IDs of the direct neighbors of unit i
#' @param  Peff 4 dimensional vector containing the estimated effects in the
#' whole population
#' @param minsize minimum number of observations for each level of the joint intervention
#' to be required in the leafs
#'
#' @return A data.frame describing the obtained Network Causal Trees.
#' Each row represents a partition (of a specific tree) with 10/14 entries.
#' Columns summary:
#' - `OF`: value of the OF in the corresponding partition,
#' - `NOBS_TR`: number of training observations in the partition,
#' - `FILTER`: values of the covariates `X` that identify the partition,
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
alleffect=function(output,tree_info,N,W,G,Y,X,Ne,Nel,p,minsize){


    if(output=="estimation"){
      data_est <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)

      NOBS_EST<-c(rep(0,nrow(tree_info)))
      EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))
      SETAU1000=SETAU1101=SETAU1110=SETAU0100=c(rep(0,nrow(tree_info)))
      tree_info<-cbind(tree_info,NOBS_EST,EFFTAU1000,SETAU1000,EFFTAU1101,SETAU1101,EFFTAU1110,SETAU1110,EFFTAU0100,SETAU0100)

      tree_info$NOBS_EST[1]<-N
      tree_info$EFFTAU1000[1]<-EffTau1000(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU1101[1]<-EffTau1101(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU1110[1]<-EffTau1110(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU0100[1]<-EffTau0100(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])

      tree_info$SETAU1000[1]<-sqrt(Vartau1000(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit],Nel=Nel))
      tree_info$SETAU1101[1]<-sqrt(Vartau1101(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit],Nel=Nel))
      tree_info$SETAU1110[1]<-sqrt(Vartau1110(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit],Nel=Nel))
      tree_info$SETAU0100[1]<-sqrt(Vartau0100(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit],Nel=Nel))

      if(nrow(tree_info)>1){
        for (j in 2:nrow(tree_info)) {

          texts=gsub(pattern="data_tree",replace="data_est",tree_info[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text=texts)))
          if(any(as.numeric(table(this_data$W,this_data$G))<minsize/2)){
            print('subpopulations not sufficiently represented in some nodes of the EST/TEST SET')
          }
          if(any(as.numeric(table(this_data$W,this_data$G))==0)){
            print('there are empty subpop in some nodes of EST/TEST SET')
          }

          Nelsub=Nel[this_data$idunit]
          tree_info$NOBS_EST[j]<-nrow(this_data)
          tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$SETAU1000[j]<- sqrt(Vartau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
          tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$SETAU1101[j]<- sqrt(Vartau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
          tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$SETAU1110[j]<- sqrt(Vartau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
          tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$SETAU0100[j]<- sqrt(Vartau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
        }}
      colnames(tree_info)<-c("OF","FILTER","TERMINAL","NOBS_TR","NOBS_EST","EFF1000_EST","SE1000_EST","EFF1101_EST","SE1101_EST","EFF1110_EST","SE1110_EST","EFF0100_EST","SE0100_EST")

    }



    if(output=="detection"){
      data_est <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)

      NOBS_EST=EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))

      tree_info<-cbind(tree_info,NOBS_EST,EFFTAU1000,EFFTAU1101,EFFTAU1110,EFFTAU0100)

      tree_info$NOBS_EST[1]<-N
      tree_info$EFFTAU1000[1]<-EffTau1000(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU1101[1]<-EffTau1101(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU1110[1]<-EffTau1110(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])
      tree_info$EFFTAU0100[1]<-EffTau0100(N=nrow(data_est),W=data_est$W,G=data_est$G,Y=data_est$Y,p=p[data_est$idunit],Ne=Ne[data_est$idunit])

      if(nrow(tree_info)>1){
        for (j in 2:nrow(tree_info)){


          texts=gsub(pattern="data_tree",replace="data_est",tree_info[j, "FILTER"])
          this_data <- subset(data_est, eval(parse(text=texts)))

          if(any(as.numeric(table(this_data$W,this_data$G))<3)){
            warning('subpopulations not sufficiently represented')
          }
          tree_info$NOBS_EST[j]<-nrow(this_data)
          tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
          tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])

        }}
      colnames(tree_info)<-c("OF","FILTER","TERMINAL","NOBS_TR","NOBS_EST","EFF1000_EST","EFF1101_EST","EFF1110_EST","EFF0100_EST")

    }
  return(tree_info)
}

