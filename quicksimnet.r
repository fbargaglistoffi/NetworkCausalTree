rm(list=ls())
library(statnet)
setwd("~/Documents/IMT/Networks and ML/Script/Github")
source("Effects.R")
source("Variances.R")
source("Probabilities.R")
source("TreeSprouting.R")
source("Others.R")
source("PlotNCT.R")
source("WholeNCT.R")

library(stringi)

library(igraph)
#set the number of nodes
N=4000
m=40

#generate covariates
x1=rbinom(N,1,0.5)
x2=rbinom(N,1,0.5)
x3=rbinom(N,1,0.3)
x4=rbinom(N,1,0.6)
X=cbind(x1,x2,x3,x4)

#Generate m networks 
gsize=N/m
library(ergm)

#if N=1500 then param ==0.035, if N=3000 then param 0.01
#if N=1500 then coefergm=c(2,-5.5), if N=3000 then coefergm=c(2,-7)
## Homophily based on x1
simnet<-genmultnet(N=N,m=m,method="ergm",varhom=x1,param = NULL, coefergm=c(2,-5.5))
#simnet<-genmultnet(N=N,m=m,method="er",varhom=x1,param = 0.01)

#group Indicator 

M=c(rep(1:m,gsize))
M=sort(M)
levels(M)<-c(1:m)
net<-graph_from_adjacency_matrix(simnet)


p=runif(m,min=0.3,max=0.3) #m dimensioned vector identifying the assignment prob. in each group

#assign individual assignment proob
prt=c()
for(i in 1:m){
  prt[which(M==i)]<-p[i]
}

#randomly assign unit to treatment arms
treat=c()
for( i in 1:N){
  treat[i]<-rbinom(1,1,prob=prt[i])
}

#generate outcome  variable
outcome <-round(rnorm(N,20,sqrt(10)),2)
#take the whole adiac_matrix
adiac_matrix<-simnet
isSymmetric(adiac_matrix)
neigh<-rowSums(adiac_matrix)

#Compute number of treated neigh and consequently Gi
num_tr_neigh <-as.vector(adiac_matrix %*% treat) 
neightreat=rep(1,N) #Gi
neightreat[which(num_tr_neigh==0)]<-0

#Pass to the standard notation
w<-treat[neigh>0]
g<-neightreat[neigh>0]
y<-outcome[neigh>0]
M<-M[neigh>0]
X<-cbind(x1,x2,x3,x4)
X<-X[neigh>0,]
N<-length(w)
adiac_matrix<-adiac_matrix[neigh>0,neigh>0]
neigh<-neigh[neigh>0]
table(w,g)

## Generate Causal Rules tau1000. Heterogeneous tau1000=i over the pop driven by x1




N=length(w)
tau=rep(2, N)
## Generate Treatment Effects
y0 = rnorm(N,sd=.1)
y1 = y0 + tau
## Generate Outcome
y = y0*(1-w) +  y1*w
hist(y)


NCT=NetworkCausalTrees(effweights<-c(1,0,0,0), 
                 A<-adiac_matrix,p<-prt,fracpredictors<-1, 
                 W<-w, Y<-y, X<-X,M<-M,  G<-g,
                 Ne<-neigh, mdisc<-25, mest<-15, 
                 minpopfrac<-1, depth<-3,
                 minsize<-5, n_trees=1, 
                 method<-"singular",output<-"detection")

title<-expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
plot.NetworkCausalTrees(NCT=NCT,output="detection",vcolor=c("seagreen4","seagreen1","lightblue1","dodgerblue2")
                        ,vsize=28,vsize2=28,coloreff ="1000",
                        vshape="circle",vframecolor = "black",
                        vlabelcolor = "black", vlabelcex =1, 
                        vlabelfont = 1,  elabelcex = 0.7, 
                        varnames=c("x1","x2","x3","x4"),
                        elabelfamily = "sans", ecolor="black",
                        ewidth=0.3, elabelcolor="black", ot=FALSE, colleg = "black", 
                        coltitleleg = "black", font.main=1,
                        cex.main = 1,adj=1,col.main = "black",
                        title=title
)



