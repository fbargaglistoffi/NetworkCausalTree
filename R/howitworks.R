rm(list=ls())

setwd("~/Documents/IMT/Networks and ML/Script/Github")
source("Effects.R")
source("Variances.R")
source("Probabilities.R")
source("TreeSprouting.R")
source("Others.R")
source("PlotNCT.R")
source("WholeNCT.R")
source("DGP.R")

library(stringi)
library(igraph)
library(statnet)
library(ergm)
library(plyr)
library(data.tree)
library(stats)
library(base)

seed=19041952



#--------------------------------------------------------------------------------------
# Singular splitting function

dgp<-NCT_data_generating_process(N=2000, m=40, p=rep(0.2,2000), het=TRUE, taui=2, K=4,
                                 method_networks="er", param_er =0.1, 
                                 var_homophily_ergm = NULL,
                                 coef_ergm = NULL)

A=dgp[[1]]
data_to_be_used=dgp[[2]]
X<-dgp[[3]]

NCT=NetworkCausalTrees(effweights<-c(1,0,0,0), 
                       A=A,p=data_to_be_used$p, fracpredictors=1, 
                       W=data_to_be_used$W, Y=data_to_be_used$Y,
                       M=as.numeric(data_to_be_used$M), G=data_to_be_used$G, X=X,
                       Ne=data_to_be_used$Ne, mdisc=25, mest=15,  
                       minpopfrac=1, depth=3,
                       minsize=5, n_trees=1, 
                       method="singular",output<-"estimation")

title<-expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
plot.NetworkCausalTrees(NCT=NCT,output="estimation",vcolor=c("seagreen4","seagreen1","lightblue1","dodgerblue2")
                        ,vsize=32,vsize2=32,coloreff ="1000",
                        vshape="circle",vframecolor = "black",
                        vlabelcolor = "black", vlabelcex =0.8, 
                        vlabelfont = 1,  elabelcex = 0.7, 
                        varnames=c("x1","x2","x3","x4"),
                        elabelfamily = "sans", ecolor="black",
                        ewidth=0.3, elabelcolor="black", ot=FALSE, colleg = "black", 
                        coltitleleg = "black", font.main=1, edge.arrow.size = 0.5,
                        cex.main = 1,adj=1,col.main = "black",
                        title=title
)


#--------------------------------------------------------------------------------------
# Composite splitting function, NCT based on all the four effects

dgp<-NCT_data_generating_process(N=2000, m=40, p=rep(0.2,2000), het=FALSE, taui= 0, K=4,
                                 method_networks="sf", param_er = NULL, 
                                 var_homophily_ergm = NULL,
                                 coef_ergm = NULL)

A=dgp[[1]]
data_to_be_used=dgp[[2]]
X<-dgp[[3]]

NCT=NetworkCausalTrees(effweights<-c(0.25,0.25,0.25,0.25), 
                       A=A,p=data_to_be_used$p, fracpredictors=1, 
                       W=data_to_be_used$W, Y=data_to_be_used$Y,
                       M=as.numeric(data_to_be_used$M), G=data_to_be_used$G, X=X,
                       Ne=data_to_be_used$Ne, mdisc=25, mest=15,  
                       minpopfrac=1, depth=3,
                       minsize=5, n_trees=1, 
                       method="composite",output<-"estimation")

title<-expression(paste("CAUSAL TREE TARGETED TO ALL THE EFFECTS"),sep="")
plot.NetworkCausalTrees(NCT=NCT,output="estimation",vcolor=c("seagreen4","seagreen1","lightblue1","dodgerblue2")
                        ,vsize=32,vsize2=32,coloreff ="1000",
                        vshape="circle",vframecolor = "black",
                        vlabelcolor = "black", vlabelcex =0.8, 
                        vlabelfont = 1,  elabelcex = 0.7, 
                        varnames=c("x1","x2","x3","x4"),
                        elabelfamily = "sans", ecolor="black",
                        ewidth=0.3, elabelcolor="black", ot=FALSE, colleg = "black", 
                        coltitleleg = "black", font.main=1, edge.arrow.size = 0.5,
                        cex.main = 1,adj=1,col.main = "black",
                        title=title
)


