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
                          K = 4,
                          m = 40, 
                          p = rep(0.2,2000), 
                          het = TRUE, 
                          taui = 2, 
                          method_networks = "er", 
                          param_er = 0.1)

result <- NetworkCausalTrees(X =  dataset[["X"]],
                             Y = dataset[["Y"]],
                             W = dataset[["W"]],
                             effweights <- c(1,0,0,0), 
                             A = dataset[["A"]],
                             G =  dataset[["G"]], 
                             M = dataset[["M"]],
                             p = dataset[["p"]], 
                             mdisc = 25, 
                             mest = 15,  
                             minpopfrac = 1, 
                             fracpredictors = 1, 
                             n_trees = 1,
                             depth = 3,
                             minsize = 5, 
                             method = "singular",
                             output = "estimation")

title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
cov_names <- colnames(dataset[["X"]])
plot(result, cov_names, title)
```

Composite splitting (NCT based on all the four effects)
```r
dataset <- data_generator(N = 2000, 
                          K = 4,
                          m = 40, 
                          p = rep(0.2,2000), 
                          het = FALSE, 
                          taui = 0, 
                          method_networks = "sf")

result <- NetworkCausalTrees(X =  dataset[["X"]],
                             Y = dataset[["Y"]],
                             W = dataset[["W"]],
                             effweights <- c(0.25,0.25,0.25,0.25), 
                             A = dataset[["A"]],
                             G =  dataset[["G"]], 
                             M = dataset[["M"]],
                             p = dataset[["p"]], 
                             mdisc = 25, 
                             mest = 15,  
                             minpopfrac = 1,
                             fracpredictors = 1, 
                             n_trees = 1, 
                             depth = 3,
                             minsize = 5, 
                             method = "composite",
                             output = "estimation")
                          
title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
cov_names <- colnames(dataset[["X"]])
plot(result, cov_names, title)
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
