#---------------------------------------------------
# 0) CLUSTERED NETWORKS GENERATION (FOR SIMULATIONS)
#--------------------------------------------------

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

#--------------------------------------
# 1) AS ESTIMATOR FUNCTIONS
#--------------------------------------


#--------------------------------------
# 1.1) PRELIMINARIES

exposure_map_AS <- function(adiac_matrix, W) {
  N <- nrow(adiac_matrix)
  G <- as.numeric(W %*% adiac_matrix)
  return(matrix(as.numeric(
    c(
      W > 0   & G > 0,
      W > 0   & G== 0,
      W == 0  & G > 0,
      W == 0  & G == 0
    )
  ),
  N, 4, dimnames = list(
    NULL, c('11', '10', '01', '00')
  )))
}

#--------------------------------------
# 1.2) EFFECTS

EffTau1000=function(N,D,Y,Pi){
  yT_HT10 <-sum(Y[D[,'10']==1] / Pi['10', D[,'10']==1])
  yT_HT00 <-sum(Y[D[,'00']==1] / Pi['00', D[,'00']==1])
  tau_HT<-(1/N)*(yT_HT10-yT_HT00)
  return(tau_HT)
}

EffTau1101=function(N,D,Y,Pi){
  yT_HT11 <-sum(Y[D[,'11']==1] / Pi['11', D[,'11']==1])
  yT_HT01 <-sum(Y[D[,'01']==1] / Pi['01', D[,'01']==1])
  tau_HT<-(1/N)*(yT_HT11-yT_HT01)
  return(tau_HT)
}

EffTau1110=function(N,D,Y,Pi){
  yT_HT11 <-sum(Y[D[,'11']==1] / Pi['11', D[,'11']==1])
  yT_HT10 <-sum(Y[D[,'10']==1] / Pi['10', D[,'10']==1])
  tau_HT<-(1/N)*(yT_HT11-yT_HT10)
  return(tau_HT)
}

EffTau0100=function(N,D,Y,Pi){
  yT_HT01 <-sum(Y[D[,'01']==1] / Pi['01', D[,'01']==1])
  yT_HT00 <-sum(Y[D[,'00']==1] / Pi['00', D[,'00']==1])
  tau_HT<-(1/N)*(yT_HT01-yT_HT00)
  return(tau_HT)
}

popeff=function(N,D,Y,Pi){
  PEffTau1000<-EffTau1000(N=N,D=D,Y=Y,Pi=Pi)
  PEffTau1101<-EffTau1101(N=N,D=D,Y=Y,Pi=Pi)  
  PEffTau1110<-EffTau1110(N=N,D=D,Y=Y,Pi=Pi) 
  PEffTau0100<-EffTau0100(N=N,D=D,Y=Y,Pi=Pi) 
  Peff<-c(PEffTau1000,PEffTau1101,PEffTau1110,PEffTau0100)
  return(Peff)
}

#--------------------------------------
# 1.3) VARIANCES

Vary <- function(w,g,D,Y,P) {
  
  nameff<-paste(w,g,sep="")
  
  #unadj part
  var <- 0
  pi <- P[[2]][[paste(nameff,nameff,sep=",")]]
  ind <- diag(pi)
  mm <- D[,nameff] %o% D[,nameff] * (pi - ind %o% ind )/pi * (Y %o% Y) / (ind %o% ind)
  mm[!is.finite(mm)] <- 0
  second_part_sum <- sum(mm)
  var <- second_part_sum
  
  #adjst
  A2 <- 0
  m <- D[,nameff]*(Y^2)/(2*ind)
  A2_part_sum <- sum(outer(m, m, FUN='+') * (pi == 0) * (!diag(length(m))) )
  A2 <- A2_part_sum
  
  #combine the two parts
  var_adjusted <- var + A2
  return(var_adjusted)
}

Covy=function(w1,g1,w2,g2,D,Y,P){
  
  #name the effects
  nameff1<-paste(w1,g1,sep="")
  nameff2<-paste(w2,g2,sep="")
  
  #first component
  cov_yT_A <- 0
  pi_k <- P[[2]][[paste(nameff1,nameff1,sep=',')]]
  ind_kk <- diag(pi_k)
  pi_l <- P[[2]][[paste(nameff2,nameff2,sep=',')]]
  ind_ll <- diag(pi_l)
  pi_k_l <-P[[3]][[ paste(nameff1 ,nameff2,sep=',') ]]
  mm <- D[,nameff1] %o% D[,nameff2] * (pi_k_l - ind_kk %o% ind_ll)/pi_k_l * (Y %o% Y) / (ind_kk %o% ind_ll)
  mm[!is.finite(mm)] <- 0
  first_part_cov <- sum(mm)
  
  #second component
  second_part_cov <- 0
  for (i in 1:nrow(D)){
    for (j in 1:nrow(D)) {
      if (pi_k_l[i,j]==0) {
        second_part_cov_i_j <- ((D[i,nameff1]*Y[i]^2/(2*ind_kk[i])) + 
                                  (D[j,nameff2]*Y[j]^2/(2*ind_ll[j])))
        second_part_cov <- second_part_cov + second_part_cov_i_j
        
      }
    }
  }
  #combine
  cov_yT_A<- first_part_cov - second_part_cov
  
  return(cov_yT_A)
}

Vartau1000=function(N,D,Y,P){
  
  vary10<-Vary(w=1,g=0,D=D,Y=Y,P=P)
  vary00<-Vary(w=0,g=0,D=D,Y=Y,P=P)
  covy10y00 <-Covy(w1=1,g1=0,w2=0,g2=0,D=D,Y=Y,P=P)
  
  var1000=abs((1/(N^2))*(vary10+vary00-2*covy10y00))
  
  return(as.numeric(var1000)) 
}

Vartau1101=function(N,D,Y,P){
  
  vary11<-Vary(w=1,g=1,D=D,Y=Y,P=P)
  vary01<-Vary(w=0,g=1,D=D,Y=Y,P=P)
  covy11y01 <-Covy(w1=1,g1=1,w2=0,g2=1,D=D,Y=Y,P=P)
  
  var1101=abs((1/(N^2))*(vary11+vary01-2*covy11y01))

  return(as.numeric(var1101)) 
}

Vartau1110=function(N,D,Y,P){
  
  vary11<-Vary(w=1,g=1,D=D,Y=Y,P=P)
  vary10<-Vary(w=1,g=0,D=D,Y=Y,P=P)
  covy11y10 <-Covy(w1=1,g1=1,w2=1,g2=0,D=D,Y=Y,P=P)
  var1110=abs((1/(N^2))*(vary11+vary10-2*covy11y10))
  
  return(as.numeric(var1110)) 
}

Vartau0100=function(N,D,Y,P){
  
  vary01<-Vary(w=0,g=1,D=D,Y=Y,P=P)
  vary00<-Vary(w=0,g=0,D=D,Y=Y,P=P)
  covy01y00 <-Covy(w1=0,g1=1,w2=0,g2=0,D=D,Y=Y,P=P)
  var0100=abs((1/(N^(2)))*(vary01+vary00-2*covy01y00))
  
  return(as.numeric(var0100)) 
}


#--------------------------------------
#  2) SPLITTING FUNCTIONS
#--------------------------------------

#--------------------------------------
#  2.0) PRELIMINARIES

get_sublist <- function(oldlist,condition) {
  newlist<-list(vector(mode = "list", length = length(names(oldlist[[1]]))) ,
                vector(mode = "list", length = length(names(oldlist[[2]]))),
                vector(mode = "list", length = length(names(oldlist[[3]]))) )
  names(newlist)=names(oldlist)
  names(newlist[[1]])<-names(oldlist[[1]])        
  names(newlist[[2]])<-names(oldlist[[2]])   
  names(newlist[[3]])<-names(oldlist[[3]])   
  names1<-names(newlist[[1]])
  names2<-names(newlist[[2]])
  names3<-names(newlist[[3]])
  for (i in names1) {
    newlist[[1]][[i]]<-oldlist[[1]][[i]][condition,]
  }
  for (i in names2) {
    newlist[[2]][[i]]<-oldlist[[2]][[i]][condition,condition]   
  }
  for (i in names3) {
    newlist[[3]][[i]]<-oldlist[[3]][[i]][condition,condition]    
  }
  return(newlist)
}

#--------------------------------------
#  2.1) COMPUTES GOF

GOF=function(method,alpha,beta,gamma,delta,N,D,Y,Pi,P,Peff,vartot,leafs){
  
  if(method=="composite"){
    ingof<-alpha*(((EffTau1000(N=N,D=D,Y=Y,Pi=Pi))^2)/(Peff[1])^2)+
           beta*(((EffTau1101(N=N,D=D,Y=Y,Pi=Pi))^2)/(Peff[2])^2)+
           gamma*(((EffTau1110(N=N,D=D,Y=Y,Pi=Pi))^2)/(Peff[3])^2)+
           delta*(((EffTau0100(N=N,D=D,Y=Y,Pi=Pi))^2)/(Peff[4])^2)  
  } 
  
  if(method=="penalized")  
  {l=leafs
  ingof<-
    alpha*(
      ((EffTau1000(N=N,D=D,Y=Y,Pi=Pi))^2) -
        2/l*sum(c(vartot,Vartau1000(N=N,D=D,Y=Y,P=P)))
    ) +
    beta*(
      ((EffTau1101(N=N,D=D,Y=Y,Pi=Pi))^2)-
        2/l*sum(c(vartot,Vartau1101(N=N,D=D,Y=Y,P=P)))
    ) *  
    gamma*(
      ((EffTau1110(N=N,D=D,Y=Y,Pi=Pi))^2)-
        2/l*sum(c(vartot,Vartau1110(N=N,D=D,Y=Y,P=P)))
    )+  
    delta*(
      ((EffTau0100(N=N,D=D,Y=Y,Pi=Pi))^2)- 
        2/l*sum(c(vartot,Vartau0100(N=N,D=D,Y=Y,P=P)))
    )
  }
  
  if (method=="singular")
  {ingof<-alpha*((EffTau1000(N=N,D=D,Y=Y,Pi=Pi))^2)+
    beta*((EffTau1101(N=N,D=D,Y=Y,Pi=Pi))^2)+
    gamma*((EffTau1110(N=N,D=D,Y=Y,Pi=Pi))^2)+
    delta*((EffTau0100(N=N,D=D,Y=Y,Pi=Pi))^2)  
  }
  return(ingof)
}

#--------------------------------------
#  2.2) SPLITS WHERE GOF IS MAXIMIZED

Gof_Split=function(method,alpha,beta,gamma,delta,N,D,Y,X,P,Pi,Peff,vartot,leafs){
  
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
      if ( all(colSums(D[x>=sp,]>0) & colSums(D[x<sp,]>0)) ) {   
        gofx[i]<- 1/2*(GOF(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta, 
                          N=length(which(x<sp)),D=D[x < sp,],Y=Y[x < sp],
                          P=get_sublist(P,which(x<sp)), Pi=Pi[,x < sp],
                          Peff=Peff,vartot=vartot,leafs=leafs) +
                      GOF(method=method,alpha=alpha,beta=beta,gamma=beta,delta=delta,
                          N=length(which(x>=sp)),D=D[x >= sp,],Y=Y[x >= sp],
                          P=get_sublist(P,which(x>=sp)), Pi=Pi[, x >= sp],
                          Peff = Peff,vartot=vartot,leafs=leafs))
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

#--------------------------------------
#  3) TRUE FUNCTIONS
#--------------------------------------

#--------------------------------------
#  3.1 BUILDS UP THE TREE DETECTING THE PARTITIONS

netctree <- function(method,alpha,beta,gamma,delta,depth,minsize,N,D,Y,X,P,Pi,Peff)
{
  
  data_tree <- data.frame(idunit=1:N,D=D,Y=Y,X=X)  
  colnames(data_tree)<-gsub(pattern = ".D", x=colnames(data_tree),replacement = "")
  do_splits <- TRUE
  
  #CREATE OUTPUT DATASET
  tree_info <- data.frame(NODE = 1, GOF=0, NOBS = nrow(data_tree), FILTER = NA, TERMINAL = "SPLIT",
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
      #SPLIT WRT GOF OVER THE SUBSET
      splitting <- Gof_Split(method=method,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                             N=nrow(this_data),Y=this_data$Y,
                             X=this_data[, grepl("X.",names(this_data))],
                             P=get_sublist(P,this_data$idunit), Pi=Pi[,this_data$idunit],
                             D=D[this_data$idunit,],
                             Peff = Peff,
                             leafs=leafs, vartot = vartot)
      
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
      
      if(method=="penalized"){
        varchild=alpha*(Vartau1000(N=nrow(this_data), Y=this_data$Y,
                                   D=D[this_data$idunit,], 
                                   P=get_sublist(P,this_data$idunit) )) +
                 beta*(Vartau1101(N=nrow(this_data), Y=this_data$Y,
                                  D=D[this_data$idunit,], 
                                  P=get_sublist(P,this_data$idunit) ))+
                 gamma*(Vartau1110(N=nrow(this_data), Y=this_data$Y,
                                   D=D[this_data$idunit,], 
                                   P=get_sublist(P,this_data$idunit) ))+
                 delta*(Vartau0100(N=nrow(this_data), Y=this_data$Y,
                                   D=D[this_data$idunit,], 
                                   P=get_sublist(P,this_data$idunit) ))
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

#--------------------------------------
#  3.2) COMPUTES THE EFFECTS


alleffect=function(output,tree_info,N,D,Y,X,Pi,P,minsize){
  
  if(output=="estimation"){
    data_est <- data.frame(idunit=1:N,D=D,Y=Y,X=X)  
    colnames(data_est)<-gsub(pattern="D.",replacement="",x=colnames(data_est))
    NOBS_EST<-c(rep(0,nrow(tree_info)))
    EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))
    SETAU1000=SETAU1101=SETAU1110=SETAU0100=c(rep(0,nrow(tree_info)))
    tree_info<-cbind(tree_info,NOBS_EST,EFFTAU1000,SETAU1000,EFFTAU1101,SETAU1101,EFFTAU1110,SETAU1110,EFFTAU0100,SETAU0100)
    
    tree_info$NOBS_EST[1]<-N
    tree_info$EFFTAU1000[1]<-EffTau1000(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU1101[1]<-EffTau1101(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU1110[1]<-EffTau1110(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU0100[1]<-EffTau0100(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])  
    
    tree_info$SETAU1000[1]<-sqrt(Vartau1000(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y, P=get_sublist(P,data_est$idunit)))
    tree_info$SETAU1101[1]<-sqrt(Vartau1101(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y, P=get_sublist(P,data_est$idunit)))
    tree_info$SETAU1110[1]<-sqrt(Vartau1110(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y, P=get_sublist(P,data_est$idunit)))
    tree_info$SETAU0100[1]<-sqrt(Vartau0100(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y, P=get_sublist(P,data_est$idunit)))
    
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

        tree_info$NOBS_EST[j]<-nrow(this_data)
        tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit] )
        tree_info$SETAU1000[j]<- sqrt(Vartau1000(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y, P=get_sublist(P,this_data$idunit)))
        tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit] )
        tree_info$SETAU1101[j]<- sqrt(Vartau1101(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y, P=get_sublist(P,this_data$idunit)))
        tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit] )
        tree_info$SETAU1110[j]<- sqrt(Vartau1110(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y, P=get_sublist(P,this_data$idunit)))
        tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit] )  
        tree_info$SETAU0100[j]<- sqrt(Vartau0100(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y, P=get_sublist(P,this_data$idunit)))
      }}
    colnames(tree_info)<-c("GOF","NOBS_TR","FILTER","N_TREE","NOBS_EST","EFF1000_EST","SE1000_EST","EFF1101_EST","SE1101_EST","EFF1110_EST","SE1110_EST","EFF0100_EST","SE0100_EST")
    
  }
  
  if(output=="detection"){
    
    data_est <- data.frame(idunit=1:N,D=D,Y=Y,X=X)  
    colnames(data_est)<-gsub(pattern="D.",replacement="",x=colnames(data_est))
    
    NOBS_EST=EFFTAU1000=EFFTAU1101=EFFTAU1110=EFFTAU0100=c(rep(0,nrow(tree_info)))
    tree_info<-cbind(tree_info,NOBS_EST,EFFTAU1000,EFFTAU1101,EFFTAU1110,EFFTAU0100)
    
    tree_info$NOBS_EST[1]<-N
    tree_info$EFFTAU1000[1]<-EffTau1000(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU1101[1]<-EffTau1101(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU1110[1]<-EffTau1110(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    tree_info$EFFTAU0100[1]<-EffTau0100(N=nrow(data_est),D=D[data_est$idunit,],Y=data_est$Y,Pi=Pi[,data_est$idunit])
    
    if(nrow(tree_info)>1){
      for (j in 2:nrow(tree_info)){
        
        
        texts=gsub(pattern="data_tree",replace="data_est",tree_info[j, "FILTER"])
        this_data <- subset(data_est, eval(parse(text=texts)))
        
        if(any(as.numeric(table(this_data$W,this_data$G))<3)){
          warning('subpopulations not sufficiently represented')  
        }
        tree_info$NOBS_EST[j]<-nrow(this_data)
        tree_info$EFFTAU1000[j]<-EffTau1000(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit])
        tree_info$EFFTAU1101[j]<-EffTau1101(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit])
        tree_info$EFFTAU1110[j]<-EffTau1110(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit])
        tree_info$EFFTAU0100[j]<-EffTau0100(N=nrow(this_data),D=D[this_data$idunit,],Y=this_data$Y,Pi=Pi[,this_data$idunit]) 
        
      }}
    colnames(tree_info)<-c("GOF","NOBS_TR","FILTER","N_TREE","NOBS_EST","EFF1000_EST","EFF1101_EST","EFF1110_EST","EFF0100_EST")
  }
  
  return(tree_info)
  
}

#-----------------------------------------------------------
#  3.3) SPROUTS THE TREE OVER A GIVRN FRAC OF UNITS AND PREDICYORS

sproutnetctree=function(method,minpopfrac,fracpredictors,sampgroup,
                        alpha,beta,gamma,delta,
                        depth,minsize,N,D,Y,X,M,Peff,P,Pi){
  
  # coerce to data.frame
  data <- data.frame(idunit=1:N,D=D,Y=Y,X=X,M=M)  
  colnames(data)<-gsub(pattern="D.",replacement="",x=colnames(data))
  #SAMPLE PREDICTORS
  samppredictors <- sort(sample(which(grepl("X.",names(data))),
                                size = ceiling(length(which(grepl("X.",names(data)))) * fracpredictors),
                                replace = FALSE))
  
  #CREATE SUBSET OF DATA ABOVE WHICH IS USED TO BUILD UP THE TREE AND BUILD THE TREE
  datasample<-data[which(M %in% sampgroup),c(1:6,samppredictors,ncol(data))]
  if(round(minpopfrac*nrow(datasample))==nrow(datasample)){
    sampunit<- nrow(datasample)
  } else {
    sampunit=sample(round(minpopfrac*nrow(datasample)):nrow(datasample),size=1)}
  
  datasample<-datasample[sample(1:nrow(datasample), size = sampunit, replace = FALSE),]  
  datasample<-datasample[order(datasample$idunit) ,]
  sampleid<-unique(datasample$idunit)
  N=length(sampleid)
  D=D[datasample$idunit,]
  Pi=Pi[,datasample$idunit]
  P=get_sublist(P,datasample$idunit)
  Y=as.numeric(datasample$Y)
  X=as.matrix(datasample[,-c(1:6,dim(datasample)[2])])
  colnames(X)=sub("X.","",colnames(X))
  
  tree_info<-netctree(method=method,alpha=alpha,beta=beta,gamma=gamma,delta=delta,N=length(sampleid),
                      depth=depth,minsize=minsize,D=D,Y=Y,X=X,
                      Pi=Pi,P=P,Peff=Peff)
  
  return(list(tree=tree_info,predictors_used=colnames(data[,samppredictors])))
}

#--------------------------------------
#  ++++ WHOLE FUNCTION ++++
#--------------------------------------


NetworkCausalTrees=function(effweights,p,A,W,Y,X,M,ngdisc,gdisc,R,
                            minpopfrac,fracpredictors,depth,minsize,n_trees,method,output){
  
   N=length(W)
   D<-exposure_map_AS(A, W=W)
   omega <- make_tr_vec_permutation( N=N, p=p,
                                     R =R, seed = 420, allow_repetitions = FALSE)
   P<- make_exposure_prob( omega, A, exposure_map_AS )
   Pi<-make_prob_exposure_cond(P)
  
  #get input weights
  alpha=effweights[1]
  beta=effweights[2]
  gamma=effweights[3]
  delta=effweights[4]
  
  
  #error messages
  if(alpha+beta+gamma+delta!=1){
    stop('weights have to sum up to one')
  }
  
  if(length(which(effweights>0))>1 & (method=="singular" | method=="penalized")){
    stop('if method is set to singular only one effect should have positive weight')
  }
  
  if(1 %in% effweights & method=="composite"){
    stop('composite gof is computed if at least two effects are investigated')
  }
  
  if( all(is.na(gdisc) & is.na(ngdisc) )){
    stop('you must specify either the groups you want to assign to the training set and estimation set, respectively
         or the number of groups that have to be sampled to be randomly assigned to the two sets')
  }
  
  data <- data.frame(idunit=1:N,D=D,Y=Y,X=X,M=M)  
  colnames(data)<-gsub(pattern="D.",replacement="",x=colnames(data))
  
  Peff<-popeff(N=N,D=D,Y=Y,Pi=Pi)
  
  m=length(unique(data$M))
  if( all(is.na(gdisc)) ) {
  sampgroup_train=sort(sample(1:m,size=ngdisc,replace=FALSE))
  }else{sampgroup_train=gdisc}
  
  trees <- plyr::raply(
    n_trees,
    sproutnetctree(method=method,sampgroup=sampgroup_train,fracpredictors=fracpredictors,minpopfrac=minpopfrac,
                   depth=depth,minsize=minsize,alpha=alpha,beta=beta,gamma=gamma,delta=delta,
                   N=N,D=D,Y=Y,X=X,M=M,Pi=Pi,P=P,Peff=Peff),
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
  sampgroup_est=sort(setdiff(1:m,sampgroup_train))
  dataest<-data[which(M %in% sampgroup_est),]  
  sampleidest<-unique(dataest$idunit)
  sampleidest<-sort(sampleidest)
  Nest=length(sampleidest)
  Dest=D[sampleidest,]
  Piest=Pi[,sampleidest]
  Pest=get_sublist(P,sampleidest)
  Yest=as.numeric(dataest$Y)
  Xest=as.matrix(dataest[,which(grepl("X.",names(dataest)))])
  colnames(Xest)=sub("X.","",colnames(Xest))

  Results_est<-alleffect(output=output,tree_info = Results,N=Nest,D=Dest,Y=Yest,X=Xest,
                         Pi=Piest,P=Pest,minsize = minsize)
  return(Results_est)
}


plot.NetworkCausalTrees=function(NCT,vcolor,vlabelcolor,
                                 ewidth,elabelcex, efamily, ecolor,elabelfamily, elabelcolor,
                                 vsize, vsize2,vshape, varnames,
                                 vlabelfont, vlabelcex,
                                 vframecolor,title,cex.main,col.main,adj){
  
  NCT$NOBS<-NCT$NOBS_EST+NCT$NOBS_TR
  #NCT$OT<-NCT$OT
  #NCT$OT_NODES_BC
  for (p in 1:length(varnames)) {
    NCT$FILTER<-gsub(pattern = paste("x",p,sep=""),replacement = varnames[p],x=as.character(NCT$FILTER) )
  }
  
  NCT$FILTER<-gsub(pattern = "data_tree",replacement ="",x=as.character(NCT$FILTER) )
  NCT$FILTER<-gsub(pattern = "[$]",replacement ="",x=as.character( NCT$FILTER ) )
  NCT$FILTER<-gsub(pattern = "_bin",replacement ="",x=as.character( NCT$FILTER ) )
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
  
  latinomefine<-rep(NULL,length(nomelati))
  edgelabel<-vector(mode = "list", length = length(nomelati))
  for (l in 1:length(nomelati)){
    if(length(grep(x=nomelati[l],pattern=">="))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[>=])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[>=])",perl = TRUE)[[1]][2:
                                                                                         (lengths(gregexpr(",", strsplit(nomelati[l],"(?=[>=])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern="<="))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[<=])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[<=])",perl = TRUE)[[1]][2:
                                                                                         (lengths(gregexpr(",", strsplit(nomelati[l],"(?=[<=])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern="<="))==0 & length(grep(x=nomelati[l],pattern="<"))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[<])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[<])",perl = TRUE)[[1]][2:
                                                                                        (lengths(gregexpr(",", strsplit(nomelati[l],"(?=[<])",perl = TRUE))) + 1)],collapse = " ")
    }
    if(length(grep(x=nomelati[l],pattern=">="))==0 & length(grep(x=nomelati[l],pattern=">"))>0){
      edgelabel[[l]][1]<-strsplit(nomelati[l],"(?=[>])",perl = TRUE)[[1]][1]
      edgelabel[[l]][2]<-paste(value=strsplit(nomelati[l],"(?=[>])",perl = TRUE)[[1]][2:
                                                                                        (lengths(gregexpr(",", strsplit(nomelati[l],"(?=[>])",perl = TRUE))) + 1)],collapse = " ")
    }
    latinomefine[l]<-paste(edgelabel[[l]][1],gsub("\\s", "", x=edgelabel[[l]][2]),sep="\n")
    
  }
  latinomefine=gsub(x=latinomefine,pattern = ">=1",replacement = "=1")
  latinomefine=gsub(x=latinomefine,pattern = "<1",replacement = "=0")
  
  E(grafo_tree)$label<-latinomefine
  eff1<-paste(round(NCT$EFF1000_EST,2),"(",round(NCT$SE1000_EST,2),")",sep="")
  eff2<-paste(round(NCT$EFF1101_EST,2),"(",round(NCT$SE1101_EST,2),")",sep="")
  eff3<-paste(round(NCT$EFF1110_EST,2),"(",round(NCT$SE1110_EST,2),")",sep="")
  eff4<-paste(round(NCT$EFF0100_EST,2),"(",round(NCT$SE0100_EST,2),")",sep="") 
  ot<-paste("(", as.numeric(round(NCT$OT,2)),")",sep="") #";",as.numeric(round(NCT$OT_NODES_BC,2)),")",sep = "")
  V(grafo_tree)$labels<-paste(eff1,eff3,NCT$NOBS,ot,sep="\n")
  NCTPLOT<-plot(grafo_tree,layout=layout_as_tree(grafo_tree), edge.label.color=elabelcolor,
                edge.width=ewidth,edge.label.cex=elabelcex, edge.label.family=elabelfamily, vertex.color=vcolor,
                vertex.label.dist = 0, vertex.label.color=vlabelcolor,
                vertex.label=V(grafo_tree)$labels,  vertex.label.font=vlabelfont, vertex.label.cex=vlabelcex,
                vertex.frame.color=vframecolor, edge.color=ecolor, vertex.frame=2,
                edge.label=E(grafo_tree)$label,vertex.shape=vshape,vertex.size=vsize,vertex.size2=vsize2)
  legend('topleft', text.col="darkblue", xjust=0, adj=0.2,
         legend=c(expression(paste(tau,"(10,00)",sep="")),expression(paste(tau,"(01,00)",sep="")),"N","OT"),title = "Legend",title.col ="red")
  title(title,cex.main=cex.main,col.main=col.main,adj=adj)
  return(NCTPLOT)
}


