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

[TODO]

## Cite

```bibtex
@article{bargagli2020heterogeneous,
  title={Heterogeneous treatment and spillover effects under clustered network interference},
  author={Bargagli-Stoffi, Falco J and Tortu, Costanza and Forastiere, Laura},
  journal={arXiv preprint arXiv:2008.00707},
  year={2020}
}
```
