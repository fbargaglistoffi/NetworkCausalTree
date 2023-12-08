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

### Example 1

Data generated using Erdos Renyi networks. 

```r
## Examples
dataset <- data_generator(N = 4000, 
                          K = 4,
                          m = 80, 
                          p = rep(0.2,4000), 
                          het = TRUE, 
                          h = 2, 
                          method_networks = "er", 
                          param_er = 0.1)
```

Singular splitting based on the main treatment effect only 

```r


result <- NetworkCausalTree(X = dataset[["X"]],
                            Y = dataset[["Y"]],
                            W = dataset[["W"]], 
                            A = dataset[["A"]],
                            M = dataset[["M"]],
                            p = dataset[["p"]], 
                            effect_weights = c(1,0,0,0),
                            ratio_dis = 0.5,
                            depth = 3,
                            minsize = 5, 
                            method = "singular",
                            output = "estimation")


title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
cov_names <- colnames(dataset[["X"]])

plot_NCT(NCT = result, 
         cov_names = cov_names,
         title = title)
```

### Example 2

Data generated using Barabasi - Albert networks. 

```r


dataset <- data_generator(N = 4000,
                          K = 4,
                          m = 80,
                          p = rep(0.2,4000),
                          het = TRUE,
                          h = 3,
                          method_networks = "sf")

```

Composite splitting (NCT based on all the four effects)

```r

result <- NetworkCausalTree(X = dataset[["X"]],
                            Y = dataset[["Y"]],
                            W = dataset[["W"]],
                            A = dataset[["A"]],
                            M = dataset[["M"]],
                            p = dataset[["p"]],
                            effect_weights =   c(0.25, 0.25, 0.25, 0.25),
                            ratio_disc = 0.5,
                            depth = 2,
                            minsize = 5,
                            method = "composite",
                            output = "detection")

title <- expression("CAUSAL TREE TARGETED TO ALL THE EFFECTS")
cov_names <- colnames(dataset[["X"]])

plot_NCT(NCT = result, 
         cov_names = cov_names,
         title = title,
         output = "detection")
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
