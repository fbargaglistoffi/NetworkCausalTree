sharedn=function(i,j,Nel){
  return(length(intersect(Nel[[i]],Nel[[j]])))
}


pi<-function(i,w,g,p,Ne){
pi<-((p[i]^w)*
    (1-p[i])^(1-w))*
    (((1-(1-p[i])^Ne[i])^g)*
    ((1-p[i])^Ne[i])^(1-g))
return(pi)
  }

pij<-function(i,j,wi,wj,gi,gj,Ne,Nel,p){
pij=((p[i]^wi)*(1-p[i])^(1-wi))*
    ((p[j]^wj)*(1-p[j])^(1-wj))*
    (((1-(1-p[i])^(Ne[i]-ifelse(Ne[i]<=Ne[j],sharedn(i,j,Nel=Nel),0)))^gi)*
    ((1-p[i])^(Ne[i])-ifelse(Ne[i]<=Ne[j],sharedn(i,j,Nel=Nel),0))^(1-gi))*
    (((1-(1-p[j])^(Ne[j]-ifelse(Ne[i]>Ne[j],sharedn(i,j,Nel=Nel),0)))^gj)*
    ((1-p[j])^(Ne[j])-ifelse(Ne[i]>Ne[j],sharedn(i,j,Nel=Nel),0))^(1-gj))
return(pij)  
}

EffTau1000=function(N,W,G,Y,p,Ne){
  tau1000=1/N*(sum(Y[W==1 & G==0]/pi(which(W==1 & G==0),1,0,p,Ne))-
               sum(Y[W==0 & G==0]/pi(which(W==0 & G==0),0,0,p,Ne)))
  return(tau1000)
}

EffTau1101=function(N,W,G,Y,p,Ne){
  tau1101=1/N*(sum(Y[W==1 & G==1]/pi(which(W==1 & G==1),1,1,p,Ne))-
          sum(Y[W==0 & G==1]/pi(which(W==0 & G==1),0,1,p,Ne)))
  return(tau1101)
}

EffTau1110=function(N,W,G,Y,p,Ne){
  tau1110=1/N*(sum(Y[W==1 & G==1]/pi(which(W==1 & G==1),1,1,p,Ne))-
          sum(Y[W==1 & G==0]/pi(which(W==1 & G==0),1,0,p,Ne)))
  return(tau1110)
}

EffTau0100=function(N,W,G,Y,p,Ne){
  tau0100=1/N*(sum(Y[W==0 & G==1]/pi(which(W==0 & G==1),0,1,p,Ne))-
          sum(Y[W==0 & G==0]/pi(which(W==0 & G==0),0,0,p,Ne)))
  return(tau0100)
}

Vartau1000=function(N,Y,W,G,p,Ne,Nel){
  
  if(length(which(W==1 & G==0))>1){
    for (i in which(W==1 & G==0)) {
    for (j in which(W==1 & G==0)) {
    if (i!=j){
      vary10=sum((1-pi(which(W==1 & G==0),1,0,p,Ne))*
             (Y[which(W==1 & G==0)]/pi(which(W==1 & G==0),1,0,p,Ne))^2) +  
             sum(((pij(i,j,1,0,1,0,Ne,Nel,p=p)-pi(i,1,0,p,Ne)*pi(j,1,0,p,Ne))/
                (pij(i,j,1,0,1,0,Ne,Nel,p=p)))*
                (Y[i]/pi(i,1,0,p,Ne))*(Y[j]/pi(j,1,0,p,Ne)) )  
    }}}}else{vary10=NA}

  if(length(which(W==0 & G==0))>1){
    for (i in which(W==0 & G==0)) {
      for (j in which(W==0 & G==0)) {
        if (i!=j){
          vary00=sum((1-pi(which(W==0 & G==0),0,0,p,Ne))*
                  (Y[which(W==0 & G==0)]/pi(which(W==0 & G==0),0,0,p,Ne))^2) +  
                  sum(((pij(i=i,j=j,0,0,0,0,Ne,Nel=Nel,p=p)-pi(i,0,0,p,Ne)*pi(j,0,0,p,Ne))/
                  (pij(i=i,j=j,0,0,0,0,Ne,Nel=Nel,p=p))) *
                  (Y[i]/pi(i,0,0,p,Ne))*(Y[j]/pi(j,0,0,p,Ne)) )  
     }}}}else{vary00=NA}
         
  if(length(which(W==1 & G==0))>1 & length(which(W==0 & G==0))>1){
  for (i in which(W==1 & G==0)) {
  for (j in which(W==0 & G==0)) {
      covy10y00=sum(1/pij(i,j,1,0,0,0,Ne,Nel,p=p)*
                (Y[i]/pi(i,1,0,p,Ne))*(Y[j]/pi(j,0,0,p,Ne)))*
                (pij(i,j,1,0,0,0,Ne,Nel,p=p)-pi(i,1,0,p,Ne)*pi(j,0,0,p,Ne))-
                (sum(((Y[which(W==1 & G==0)])^2)/(2*pi(which(W==1 & G==0),1,0,p,Ne))) +
                 sum(((Y[which(W==0 & G==0)])^2)/(2*pi(which(W==0 & G==0),0,0,p,Ne))))  
  }}}else{covy10y00=NA}
  
  if(any(is.na(c(vary10,vary00)))){
    var1000=NA
  }else{  
var1000=(1/N^{2})*(vary10+vary00-2*covy10y00)} 
return(var1000) 
}

Vartau1101=function(N,Y,W,G,p,Ne,Nel){
  
  if(length(which(W==1 & G==1))>1){
    for (i in which(W==1 & G==1)) {
      for (j in which(W==1 & G==1)) {
        if (i!=j){
          vary11=sum((1-pi(which(W==1 & G==1),1,1,p,Ne))*
                       (Y[which(W==1 & G==1)]/pi(which(W==1 & G==1),1,1,p,Ne))^2) +  
            sum( ((pij(i,j,1,1,1,1,Ne,Nel,p=p)-pi(i,1,1,p,Ne)*pi(j,1,1,p,Ne))/
                   (pij(i,j,1,1,1,1,Ne,Nel,p=p)))*
                   (Y[i]/pi(i,1,1,p,Ne))*(Y[j]/pi(j,1,1,p,Ne)) )  
        }}}}else{vary11=NA}
  
  if(length(which(W==0 & G==1))>1){
  for (i in which(W==0 & G==1)) {
    for (j in which(W==0 & G==1)) {
      if (i!=j){
        vary01=sum((1-pi(which(W==0 & G==1),0,1,p,Ne))*
                  (Y[which(W==0 & G==1)]/pi(which(W==0 & G==1),0,1,p,Ne))^2) +  
               sum( ((pij(i,j,0,1,0,1,Ne,Nel,p)-pi(i,0,1,p,Ne)*pi(j,0,1,p,Ne))/
                    (pij(i,j,0,1,0,1,Ne,Nel,p))) *
                    (Y[i]/pi(i,0,1,p,Ne))*(Y[j]/pi(j,0,1,p,Ne)) )  
      }}}}else{vary01=NA}
 
  if(length(which(W==1 & G==1))>1 & length(which(W==0 & G==1))>1){
  for (i in which(W==1 & G==1)) {
    for (j in which(W==0 & G==1)) {
      covy11y01=sum(1/pij(i,j,1,1,0,1,Ne,Nel,p=p)*
                (Y[i]/pi(i,1,1,p,Ne))*(Y[j]/pi(j,0,1,p,Ne)))*
                (pij(i,j,1,1,0,1,Ne,Nel,p=p)-pi(i,1,1,p,Ne)*pi(j,0,1,p,Ne)) -
                (sum(((Y[which(W==1 & G==1)])^2)/(2*pi(which(W==1 & G==1),1,1,p,Ne))) +
                sum(((Y[which(W==0 & G==1)])^2)/(2*pi(which(W==0 & G==1),0,1,p,Ne))))  
    }}}else{covy11y01=NA}
  
  if(any(is.na(c(vary11,vary01)))){
    var1101=NA
  }else{
  var1101=(1/N^{2})*(vary11+vary01-2*covy11y01)}
  return(var1101) 
}

Vartau1110=function(N,Y,W,G,p,Ne,Nel){
  
  if(length(which(W==1 & G==1))>1){
  for (i in which(W==1 & G==1)) {
    for (j in which(W==1 & G==1)) {
      if (i!=j){
        vary11=sum((1-pi(which(W==1 & G==1),1,1,p,Ne))*
               (Y[which(W==1 & G==1)]/pi(which(W==1 & G==1),1,1,p,Ne))^2) +  
          sum( ((pij(i,j,1,1,1,1,Ne,Nel,p=p)-pi(i,1,1,p,Ne)*pi(j,1,1,p,Ne))/
                 (pij(i,j,1,1,1,1,Ne,Nel,p=p)))*
                 (Y[i]/pi(i,1,1,p,Ne))*(Y[j]/pi(j,1,1,p,Ne)) )  
      }}}}else{vary11=NA}
  
  if(length(which(W==1 & G==0))>1){
  for (i in which(W==1 & G==0)) {
    for (j in which(W==1 & G==0)) {
      if (i!=j){
        vary10=sum((1-pi(which(W==1 & G==0),1,0,p,Ne))*
                     (Y[which(W==1 & G==0)]/pi(which(W==1 & G==0),1,0,p,Ne))^2) +  
          sum( ((pij(i,j,1,0,1,0,Ne,Nel,p)-pi(i,1,0,p,Ne)*pi(j,1,0,p,Ne))/
                 (pij(i,j,1,0,1,0,Ne,Nel,p)))*
                 (Y[i]/pi(i,1,0,p,Ne))*(Y[j]/pi(j,1,0,p,Ne)) )  
      }}}}else{vary10=NA}
  
  if(length(which(W==1 & G==1))>1 & length(which(W==1 & G==0))>1){
  for (i in which(W==1 & G==1)) {
    for (j in which(W==1 & G==0)) {
      covy11y10=sum(1/pij(i,j,1,1,1,0,Ne,Nel,p=p)*
                (Y[i]/pi(i,1,1,p,Ne))*(Y[j]/pi(j,1,0,p,Ne)))*
                (pij(i,j,1,1,1,0,Ne,Nel,p=p)-pi(i,1,1,p,Ne)*pi(j,1,0,p,Ne))-
                (sum(((Y[which(W==1 & G==1)])^2)/(2*pi(which(W==1 & G==1),1,1,p,Ne))) +
                sum(((Y[which(W==1 & G==0)])^2)/(2*pi(which(W==1 & G==0),1,0,p,Ne))))  
    }}}else{covy11y10=NA}
  
  if(any(is.na(c(vary11,vary10)))){
  var1110=NA
  }else{
  var1110=(1/N^{2})*(vary11+vary10-2*covy11y10)}  
  return(var1110) 
}

Vartau0100=function(N,Y,W,G,p,Ne,Nel){
  
  if(length(which(W==0 & G==1))>1){
  for (i in which(W==0 & G==1)) {
    for (j in which(W==0 & G==1)) {
      if (i!=j){
        vary01=sum((1-pi(which(W==0 & G==1),0,1,p,Ne))*
                     (Y[which(W==0 & G==1)]/pi(which(W==0 & G==1),0,1,p,Ne))^2) +  
          sum( ((pij(i,j,0,1,0,1,Ne,Nel,p=p)-pi(i,0,1,p,Ne)*pi(j,0,1,p,Ne))/
                (pij(i,j,0,1,0,1,Ne,Nel,p=p)))*
                (Y[i]/pi(i,0,1,p,Ne))*(Y[j]/pi(j,0,1,p,Ne)) )  
      }}}}else{vary01=NA}
  
  if(length(which(W==0 & G==0))>1){
  for (i in which(W==0 & G==0)) {
    for (j in which(W==0 & G==0)) {
      if (i!=j){
        vary00=sum((1-pi(which(W==0 & G==0),0,0,p,Ne))*
                     (Y[which(W==0 & G==0)]/pi(which(W==0 & G==0),0,0,p,Ne))^2) +  
          sum( ((pij(i=i,j=j,0,0,0,0,Ne,Nel=Nel,p=p)-pi(i,0,0,p,Ne)*pi(j,0,0,p,Ne))/
               (pij(i=i,j=j,0,0,0,0,Ne,Nel=Nel,p=p))) *
               (Y[i]/pi(i,0,0,p,Ne))*(Y[j]/pi(j,0,0,p,Ne)) )  
      }}}}else{vary00=NA}
  
  if(length(which(W==0 & G==1))>1 & length(which(W==0 & G==0))>1){
   for (i in which(W==0 & G==1)) {
    for (j in which(W==0 & G==0)) {
      covy01y00=sum(1/pij(i,j,0,1,0,0,Ne,Nel,p=p)*
                (Y[i]/pi(i,0,1,p,Ne))*(Y[j]/pi(j,0,0,p,Ne)))*
                (pij(i=i,j=j,0,1,0,0,Ne=Ne,Nel=Nel,p=p)-pi(i,0,1,p,Ne)*pi(j,0,0,p,Ne)) -
                (sum(((Y[which(W==0 & G==1)])^2)/(2*pi(which(W==0 & G==1),0,1,p,Ne))) +
                 sum(((Y[which(W==0 & G==0)])^2)/(2*pi(which(W==0 & G==0),0,0,p,Ne))))  
    }}}else{covy01y00=NA}
  
  if(any(is.na(c(vary01,vary00)))){
  var0100=NA
  }else{
  var0100=(1/N^{2})*(vary01+vary00-2*covy01y00)}
  return(var0100) 
}

popeff=function(N,W,G,Y,p,Ne){
  PEffTau1000<-EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  PEffTau1101<-EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)  
  PEffTau1110<-EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne) 
  PEffTau0100<-EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne) 
  Peff<-c(PEffTau1000,PEffTau1101,PEffTau1110,PEffTau0100)
  return(Peff)
}

######---->SPLITTING FUNCTIONS 

######Function thst computes the composite GOF (objective function)

GOF=function(method,alpha,beta,gamma,delta,N,W,G,Y,p,Ne,Peff){
  
  if(method=="composite"){
    ingof<-alpha*(((EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[1])^2)+
            beta*(((EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[2])^2)+
            gamma*(((EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[3])^2)+
            delta*(((EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)/(Peff[4])^2)  
  } 
  if (method=="singular")
  {ingof<-alpha*((EffTau1000(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
          beta*((EffTau1101(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
          gamma*((EffTau1110(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)+
          delta*((EffTau0100(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne))^2)  
            }
  return(ingof)
}

Gof_Split=function(method,alpha,beta,gamma,delta,N,W,G,Y,X,p,Ne,Peff){
  
  #initialize
  gof <- c()
  name <-c()
  values<-c()
  #check all the predictors
  for(j in 1:dim(X)[2]){
    x=X[,j]
    splits <- sort(unique(x))
    valuesx<-c()
    gofx <- c()
    namesx<-c()
    for(i in seq_along(splits)){
      sp <- splits[i]
      if (all(as.numeric(table(W[x>=sp],G[x>=sp]))>2) & all(as.numeric(table(W[x<sp],G[x<sp]))>2)) {   
        gofx[i]<- 1/2*(GOF(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta,N=length(which(x<sp)),W=W[x < sp],G=G[x < sp],Y=Y[x < sp],Ne=Ne[x < sp],p=p[x < sp],Peff=Peff) +
                       GOF(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta,N=length(which(x>=sp)),W=W[x >= sp],G=G[x >= sp],Y=Y[x >= sp],Ne=Ne[x >= sp],p=p[x >= sp], Peff = Peff))
      } else {gofx[i]<-0}
    }
    namex=rep(colnames(X)[j],length(unique(x)))
    valuesx=c(sort(unique(x)))
    
    #append all the computed values of the GOF, all the values that have defined the split and the name of the variable that has been used.
    gof=c(gof,gofx)
    name=c(name,namex)
    values=c(values,valuesx)
  }

  if(all(is.na(gof))){
    warning('gof returns all NANs')
    gofres<-NA
    splitres<-NA
    varres<-NA
  }else{
    gofres<-max(na.omit(gof))
    splitres<-values[which.max(gof)]
    varres<-name[which.max(gof)]
  }
  

  #find the minimum value of the gof and returns the variable and the value where the minimum is reached.  
  return(c(gof = gofres , split = splitres , var=varres))

  }

netctree <- function(method,alpha,beta,gamma,delta,depth,minsize,N,W,G,Y,X,p,Ne,Peff)
{
  
  data_tree <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)  
  do_splits <- TRUE
  
  #CREATE OUTPUT DATASET
  tree_info <- data.frame(NODE = 1, GOF=0, NOBS = nrow(data_tree), FILTER = NA, TERMINAL = "SPLIT",
                          stringsAsFactors = FALSE)
  
  while(do_splits){
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")
    
    for (j in to_calculate) {
      
      if (!is.na(tree_info[j, "FILTER"])){
        texts=tree_info[j, "FILTER"]
        this_data <- subset(data_tree, eval(parse(text=texts)))
      } else {
        this_data<- data_tree}
      
      #SPLIT WRT GOF OVER THE SUBSET
      splitting <- Gof_Split(method=method,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                             N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,
                             X=this_data[, grepl("X.",names(this_data))],
                             Ne=Ne[this_data$idunit],p=p[this_data$idunit],
                             Peff = Peff)
      
      if (any(is.na(splitting))) {
        split_here <- rep(FALSE, 2)
        print('splits has stopped couse gof is all NA')
      }else{
      
      #GET THE MAX GOF
      maxgof <- as.numeric(splitting[1])
      mn <- max(tree_info$NODE)
      
      
      #PASTE FILTER RULES
      tmp_filter <- c(paste("data_tree$",splitting[3], ">=", 
                            as.numeric(splitting[2]),sep=""),
                      paste("data_tree$",splitting[3], "<", 
                            as.numeric(splitting[2]),sep=""))
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
      
      #STOP IF THE THREE IS OVER THE DEPTH
      depth_tree<-as.numeric(stri_count_regex(tree_info[j, "FILTER"], "X."))
      if (depth_tree>=depth & !is.na(depth_tree)){
        split_here <- rep(FALSE, 2)
        print('split has stopped for reached depth')
      }
      
      #STOP IF INSUFFICIENT MINSIZE
      if (any(as.numeric(table(this_data$W,this_data$G))< minsize)) {
        split_here <- rep(FALSE, 2)
        print('split has stopped for insufficient minsize')
      }
      
      #CREATE CHILDREN DATASET
      children <- data.frame(NODE = c(mn+1, mn+2),
                             GOF=c(rep(maxgof,2)),
                             NOBS = tmp_nobs,
                             FILTER = tmp_filter,
                             TERMINAL = rep("SPLIT", 2),
                             row.names = NULL)[split_here,]
      
      #ASSIGN PARENTS OR LEAF IF SPLIT_HERE IS ACTIVE OR NOT 
      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")
      tree_info <- rbind(tree_info, children) 
      #STOP SPLITTING IF THERE AREN'T NODES WITH "SPLIT"
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    } 
    
  } 
  
  return(tree = tree_info)
}


alleffect=function(output,tree_info,N,W,G,Y,X,Ne,Nel,p,minsize){

  if(output=="estimation"){
  data_est <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)  
  
  EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))
  SETAU1000=SETAU1101=SETAU1110=SETAU0100=c(rep(0,nrow(tree_info)))
  tree_info<-cbind(tree_info,EFFTAU1000,SETAU1000,EFFTAU1101,SETAU1101,EFFTAU1110,SETAU1110,EFFTAU0100,SETAU0100)

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
      warning('subpopulations not sufficiently represented')  
    }
    
    Nelsub=Nel[this_data$idunit]
   
    tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
    tree_info$SETAU1000[j]<- sqrt(Vartau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
    tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
    tree_info$SETAU1101[j]<- sqrt(Vartau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
    tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
    tree_info$SETAU1110[j]<- sqrt(Vartau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
    tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])  
    tree_info$SETAU0100[j]<- sqrt(Vartau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit],Nel=Nelsub))
  }}
  colnames(tree_info)<-c("GOF","NOBS","FILTER","NUMTREE","EFF1000_EST","SE1000_EST","EFF1101_EST","SE1101_EST","EFF1110_EST","SE1110_EST","EFF0100_EST","SE0100_EST")
  
  }

  if(output=="detection"){
    data_est <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X)  
    
    EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))
    tree_info<-cbind(tree_info,EFFTAU1000,EFFTAU1101,EFFTAU1110,EFFTAU0100)
    
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
        
        tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
        tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
        tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit])
        tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),W=this_data$W,G=this_data$G,Y=this_data$Y,p=p[this_data$idunit],Ne=Ne[this_data$idunit]) 
        
      }}
    colnames(tree_info)<-c("GOF","NOBS","FILTER","NUMTREE","EFF1000_EST","EFF1101_EST","EFF1110_EST","EFF0100_EST")
     }
  
  return(tree_info)
  
}

sproutnetctree=function(method,minpopfrac,fracpredictors,sampgroup,m,alpha,beta,gamma,delta,depth,minsize,N,W,G,Y,X,M,p,Ne,Peff){
  
  # coerce to data.frame
  data <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X,M=M)  
  
  #SAMPLE PREDICTORS
  samppredictors <- sort(sample(which(grepl("X.",names(data))),
                                size = ceiling(length(which(grepl("X.",names(data)))) * fracpredictors),
                                replace = FALSE))
  
  #CREATE SUBSET OF DATA ABOVE WHICH IS USED TO BUILD UP THE TREE AND BUILD THE TREE
  datasample<-data[which(M %in% sampgroup),c(1:4,samppredictors,ncol(data))]
  if(round(minpopfrac*nrow(datasample))==nrow(datasample)){
    sampunit<- nrow(datasample)
  } else {
    sampunit=sample(round(minpopfrac*nrow(datasample)):nrow(datasample),size=1)}
  
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
                      p=p[sampleid],Ne=Ne[sampleid],Peff=Peff)
  
  return(list(tree=tree_info,predictors_used=colnames(data[,samppredictors])))
}


NetworkCausalTrees=function(effweights,A,p,fracpredictors,W,Y,X,M,G,Ne,mdisc,mest,minpopfrac,depth,minsize,n_trees,method,output){
  
  N<-length(W)
  m<-length(unique(M))
  
  #get input weights
  alpha=effweights[1]
  beta=effweights[2]
  gamma=effweights[3]
  delta=effweights[4]
  
  
  #error messages
  if(alpha+beta+gamma+delta!=1){
    stop('weights have to sum up to one')
  }
  
  if(length(which(effweights>0))>1 & method=="singular"){
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
  
  if(output=="estimation"){
    if(is.null(G)){
      nt <-as.vector(A %*% W) 
      G=rep(1,N) 
      G[nt==0]<-0}
     Nel<-vector(mode = "list", length = N)
    for( i in  1:N){
      Nel[[i]]<-which(A[i,]>0)  
    }
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
                   N=N,W=W,G=G,Y=Y,X=X,M=M,Ne=Ne,p=p,Peff=Peff),
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
  
  Results<-data.frame(GOF=Forest$GOF,NOBS=Forest$NOBS ,  FILTER =  c("NA",as.vector(na.omit(unique(Forest$FILTER)))), NUMTREE=NA,
                      stringsAsFactors = FALSE)
  
  for(j in as.vector(na.omit(unique(Forest$FILTER)))){
    Results$NUMTREE[Results$FILTER==j]<-length(unique(Forest$IDTREE[which(Forest$FILTER==j)]))
    Results$NUMTREE[Results$FILTER=="NA"]=length(unique(Forest$IDTREE[is.na(Forest$FILTER)]))
  }  

  if(nrow(Results)==1){
    warning('No split has been made')
  } 
  # coerce to data.frame
  #estimation set
  sampgroup_est=sample(setdiff(1:m,sampgroup_train),size=mest,replace=FALSE) 
  dataest<-data[which(M %in% sampgroup_est),]  
  sampleidest<-unique(dataest$idunit)
  sampleidest<-sort(sampleidest)
  Nest=length(sampleidest)
  West=as.numeric(dataest$W)
  Gest=as.numeric(dataest$G)
  Yest=as.numeric(dataest$Y)
  Xest=as.matrix(dataest[,which(grepl("X.",names(dataest)))])
  colnames(Xest)=sub("X.","",colnames(Xest))
  if(output=="estimation"){
  Nelest=Nel[sampleidest]}
  Neighest=Ne[sampleidest]
  pest=p[sampleidest]
  
  Results_est<-alleffect(output=output,tree_info = Results,N=Nest,W=West,G=Gest,Y=Yest,X=Xest,Ne=Neighest,p=p[sampleidest],Nel=Nelest,minsize=minsize)
   return(Results_est)
}
 




SimNetworkCausalTrees=function(effweights,A,p,fracpredictors,W,Y,X,M,G,Ne,mdisc,mest,minpopfrac,depth,minsize,n_trees,method,output){
  
  N<-length(W)
  m<-length(unique(M))
  
  #get input weights
  alpha=effweights[1]
  beta=effweights[2]
  gamma=effweights[3]
  delta=effweights[4]
  
  
  #error messages
  if(alpha+beta+gamma+delta!=1){
    stop('weights have to sum up to one')
  }
  
  if(length(which(effweights>0))>1 & method=="singular"){
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
    stop('if output is set to estimation you have to specify the adiacency matrix')
  }
  
  if(is.null(Ne)){
    Ne<-rowSums(A)
  }
  
  if(output=="estimation"){
    if(is.null(G)){
      nt <-as.vector(A %*% W) 
      G=rep(1,N) 
      G[nt==0]<-0}
      Nel<-vector(mode = "list", length = N)
      for( i in  1:N){
        Nel[[i]]<-which(A[i,]>0)  
      }
  }
  

  data <- data.frame(idunit=1:N,W=W,G=G,Y=Y,X=X,M=M)  
  
  
  Peff<-popeff(N=N,W=W,G=G,Y=Y,p=p,Ne=Ne)
  
  sampgroup_train=sample(1:m,size=mdisc,replace=FALSE) 
  
  trees <- plyr::raply(
    n_trees,
    sproutnetctree(method=method,sampgroup=sampgroup_train,fracpredictors=fracpredictors,m=m,minpopfrac=minpopfrac,
                   depth=depth,minsize=minsize,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                   N=N,W=W,G=G,Y=Y,X=X,M=M,Ne=Ne,p=p,Peff=Peff),
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
  
  Results<-data.frame(GOF=Forest$GOF,NOBS=Forest$NOBS ,  FILTER =  c("NA",as.vector(na.omit(unique(Forest$FILTER)))), NUMTREE=NA,
                      stringsAsFactors = FALSE)
  
  for(j in as.vector(na.omit(unique(Forest$FILTER)))){
    Results$NUMTREE[Results$FILTER==j]<-length(unique(Forest$IDTREE[which(Forest$FILTER==j)]))
    Results$NUMTREE[Results$FILTER=="NA"]=length(unique(Forest$IDTREE[is.na(Forest$FILTER)]))
  }  
  
  if(nrow(Results)==1){
    warning('No split has been made')
  }
  
  # coerce to data.frame
  #estimation set
  sampgroup_est=sample(setdiff(1:m,sampgroup_train),size=mest,replace=FALSE) 
  dataest<-data[which(M %in% sampgroup_est),]  
  sampleidest<-unique(dataest$idunit)
  sampleidest<-sort(sampleidest)
  Nest=length(sampleidest)
  West=as.numeric(dataest$W)
  Gest=as.numeric(dataest$G)
  Yest=as.numeric(dataest$Y)
  Xest=as.matrix(dataest[,which(grepl("X.",names(dataest)))])
  colnames(Xest)=sub("X.","",colnames(Xest))
  
  if(output=="estimation"){
  Nelest=Nel[sampleidest]}
  Neighest=Ne[sampleidest]
  pest=p[sampleidest]
  
  Results_est<-alleffect(output=output,tree_info = Results,N=Nest,W=West,G=Gest,Y=Yest,X=Xest,Ne=Neighest,p=p[sampleidest],Nel=Nelest,minsize=minsize)

  #testing set
  sampgroup_test=setdiff(1:m,c(sampgroup_train,sampgroup_est))
  datatest<-data[which(M %in% sampgroup_test),]  
  sampleidtest<-unique(datatest$idunit)
  sampleidtest<-sort(sampleidtest)
  Ntest=length(sampleidtest)
  Wtest=as.numeric(datatest$W)
  Gtest=as.numeric(datatest$G)
  Ytest=as.numeric(datatest$Y)
  Xtest=as.matrix(datatest[,which(grepl("X.",names(datatest)))])
  colnames(Xtest)=sub("X.","",colnames(Xtest))
  if(output=="estimation"){
  Neltest=Nel[sampleidtest]}
  Neightest=Ne[sampleidtest]
  ptest=p[sampleidtest]
  

  Results_test<-alleffect(output=output,tree_info = Results,N=Ntest,W=Wtest,G=Gtest,Y=Ytest,X=Xtest,Ne=Neightest,p=ptest,Nel=Neltest,minsize=minsize)
  
  colnames(Results_test)<-gsub(colnames(Results_test),pattern = "_EST",replacement = "_TEST")
  
  Results<-cbind(Results,Results_est[,-c(1:4)],Results_test[,-c(1:4)])
  return(Results) 

}

plot.NetworkCausalTrees=function(NCT,vcolor,vlabelcolor,
                                 ewidth,elabelcex, efamily, ecolor,elabelfamily, elabelcolor,
                                 vsize, vsize2,vshape,
                                 vlabelfont, vlabelcex,
                                 vframecolor,title,cex.main,col.main,adj){
  
  NCT$FILTER<-gsub(pattern = "data_tree",replacement ="",x=as.character(NCT$FILTER) )
  NCT$FILTER<-gsub(pattern = "[$]",replacement ="",x=as.character( NCT$FILTER ) )
  NCT$FILTER[which(NCT$FILTER!="NA")]<-paste0(" & ", NCT$FILTER[which(NCT$FILTER!="NA")])
  NCT$FILTER[which(NCT$FILTER!="NA")]<-paste0(" NA ",NCT$FILTER[which(NCT$FILTER!="NA")])
  
  tree_data<-as.Node(NCT, mode = "table",
                     pathName = "FILTER", pathDelimiter = " & ", colLevels = NULL,
                     na.rm = TRUE)
  lati_tree<-ToDataFrameNetwork(tree_data)
  lati_tree<-as.matrix(lati_tree)
  nomelati<-NULL
  for(i in 1:nrow(lati_tree)){
    nomelati_i<-tail(strsplit(lati_tree[,2], "/")[[i]],n=1)
    nomelati_i<-gsub(pattern = "X.",replacement="",x=nomelati_i)
    nomelati=c(nomelati,nomelati_i)
  } 
  grafo_tree<-graph_from_edgelist(lati_tree, directed = TRUE)
  V(grafo_tree)$TAU1000<-NCT$EFF1000_EST
  V(grafo_tree)$SE1000<-NCT$SE1000_EST
  V(grafo_tree)$TAU1101<-NCT$EFF1101_EST
  V(grafo_tree)$SE1101<-NCT$SE1101_EST
  V(grafo_tree)$TAU1110<-NCT$EFF1110_EST
  V(grafo_tree)$SE1110<-NCT$SE1110_EST
  V(grafo_tree)$TAU0100<-NCT$EFF0100_EST
  V(grafo_tree)$SE0100<-NCT$SE0100_EST
  E(grafo_tree)$label<-nomelati
  
  
  eff1<-paste(round(NCT$EFF1000_EST,2),"(",round(NCT$SE1000_EST,2),")",sep="")
  eff2<-paste(round(NCT$EFF1101_EST,2),"(",round(NCT$SE1101_EST,2),")",sep="")
  eff3<-paste(round(NCT$EFF1110_EST,2),"(",round(NCT$SE1110_EST,2),")",sep="")
  eff4<-paste(round(NCT$EFF0100_EST,2),"(",round(NCT$SE0100_EST,2),")",sep="")  
  V(grafo_tree)$labels<-paste(eff1,eff2,eff3,eff4,sep="\n")
  NCTPLOT<-plot(grafo_tree,layout=layout_as_tree(grafo_tree), edge.label.color=elabelcolor,
                edge.width=ewidth,edge.label.cex=elabelcex, edge.label.family=elabelfamily, vertex.color=vcolor,
                vertex.label.dist = 0, vertex.label.color=vlabelcolor,
                vertex.label=V(grafo_tree)$labels,  vertex.label.font=vlabelfont, vertex.label.cex=vlabelcex,
                vertex.frame.color=vframecolor, edge.color=ecolor,
                edge.label=E(grafo_tree)$label,vertex.shape=vshape,vertex.size=vsize,vertex.size2=vsize2)
  legend('topleft', text.col="darkblue", xjust=0, adj=0.3,
         legend=c("tau(10,00)","tau(11,01)","eta(11,10)","eta(01,00)"),title = "Effects",title.col ="red")
  title(title,cex.main=cex.main,col.main=col.main,adj=adj)
  return(NCTPLOT)
}

