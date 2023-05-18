# Network Causal Tree

The `NetworkCausalTree` package introduces a machine learning method that uses tree-based algorithms and an Horvitz-Thompson estimator to assess the heterogeneity of treatment and spillover effects in clustered network interference. Causal inference studies typically assume no interference between individuals, but in real-world scenarios where individuals are interconnected through social, physical, or virtual ties, the effect of a treatment can spill over to other connected individuals in the network. To avoid biased estimates of treatment effects, interference should be considered. Understanding the heterogeneity of treatment and spillover effects can help policy-makers scale up interventions, target strategies more effectively, and generalize treatment spillover effects to other populations.

## Getting Started

Installing the latest developing version: 

```r
library(devtools)
install_github("fbargaglistoffi/NetworkCausalTree", ref="master")
```

Import:

```r
library("NetworkCausalTree")
```

## Examples

Singular splitting
```r
dataset <- data_generator(N = 2000, 
                          m = 40, 
                          p = rep(0.2,2000), 
                          het = TRUE, 
                          taui = 2, 
                          K = 4,
                          method_networks = "er", 
                          param_er = 0.1)

NCT <- NetworkCausalTrees(effweights <- c(1,0,0,0), 
                          A = dataset[["adj_matrix"]],
                          p = dataset[["p"]], 
                          fracpredictors = 1, 
                          W = dataset[["W"]],
                          Y = dataset[["Y"]],
                          M = dataset[["M"]],
                          G =  dataset[["G"]], 
                          X =  dataset[["X"]],
                          Ne =  dataset[["Ne"]],
                          mdisc = 25, 
                          mest = 15,  
                          minpopfrac = 1, 
                          depth = 3,
                          minsize = 5, 
                          n_trees = 1, 
                          method = "singular",
                          output = "estimation")

title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
plot.NetworkCausalTrees(NCT = NCT,
                        output = "estimation",
                        vcolor = c("seagreen4","seagreen1","lightblue1","dodgerblue2"),
                        vsize = 32,
                        vsize2 = 32,
                        coloreff ="1000",
                        vshape = "circle",
                        vframecolor = "black",
                        vlabelcolor = "black", 
                        vlabelcex = 0.8, 
                        vlabelfont = 1,  
                        elabelcex = 0.7, 
                        varnames = c("x1","x2","x3","x4"),
                        elabelfamily = "sans", 
                        ecolor = "black",
                        ewidth = 0.3, 
                        elabelcolor = "black", 
                        ot = FALSE, 
                        colleg = "black", 
                        coltitleleg = "black", 
                        font.main = 1, 
                        edge.arrow.size = 0.5,
                        cex.main = 1,
                        adj = 1,
                        col.main = "black",
                        title = title)
```

Composite splitting (NCT based on all the four effects)
```r
dataset <- data_generator(N = 2000, 
                          m = 40, 
                          p = rep(0.2,2000), 
                          het = FALSE, 
                          taui = 0, 
                          K = 4,
                          method_networks = "sf")

NCT <- NetworkCausalTrees(effweights <- c(0.25,0.25,0.25,0.25), 
                          A = dataset[["adj_matrix"]],
                          p = dataset[["p"]], 
                          fracpredictors = 1, 
                          W = dataset[["W"]],
                          Y = dataset[["Y"]],
                          M = dataset[["M"]],
                          G =  dataset[["G"]], 
                          X =  dataset[["X"]],
                          Ne =  dataset[["Ne"]],
                          mdisc = 25, 
                          mest = 15,  
                          minpopfrac = 1, 
                          depth = 3,
                          minsize = 5, 
                          n_trees = 1, 
                          method = "composite",
                          output = "estimation")
                          
title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
plot.NetworkCausalTrees(NCT = NCT,
                        output = "estimation",
                        vcolor = c("seagreen4","seagreen1","lightblue1","dodgerblue2"),
                        vsize = 32,
                        vsize2 = 32,
                        coloreff ="1000",
                        vshape = "circle",
                        vframecolor = "black",
                        vlabelcolor = "black", 
                        vlabelcex = 0.8, 
                        vlabelfont = 1,  
                        elabelcex = 0.7, 
                        varnames = c("x1","x2","x3","x4"),
                        elabelfamily = "sans", 
                        ecolor = "black",
                        ewidth = 0.3, 
                        elabelcolor = "black", 
                        ot = FALSE, 
                        colleg = "black", 
                        coltitleleg = "black", 
                        font.main = 1, 
                        edge.arrow.size = 0.5,
                        cex.main = 1,
                        adj = 1,
                        col.main = "black",
                        title = title)
```


## Cite

```bibtex
@article{bargagli2020heterogeneous,
  title={Heterogeneous treatment and spillover effects under clustered network interference},
  author={Bargagli-Stoffi, Falco J and Tortu, Costanza and Forastiere, Laura},
  journal={arXiv preprint arXiv:2008.00707},
  year={2020}
}
```
