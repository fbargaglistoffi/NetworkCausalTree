# Network Causal Tree algorithm

The bulk of causal inference studies rules out the presence of interference between units. However, in many real-world scenarios units are interconnected by social, physical or virtual ties and the effect of a treatment can spill from one unit to other connected individuals in the network. In these settings, interference should be taken into account to avoid biased estimates of the treatment effect, but it can also be leveraged to save resources and provide the intervention to a lower percentage of the population where the treatment is more effective and where the effect can spill over to other susceptible individuals. In fact, different people might respond differently not only to the treatment received but also to the treatment received by their network contacts. Understanding the heterogeneity of treatment and spillover effects can help policy-makers in the scale-up phase of the intervention, it can guide the design of targeting strategies with the ultimate goal of making the interventions more cost-effective, and it might even allow generalizing the level of treatment spillover effects in other populations. In this repository, we develop a machine learning method that makes use of tree-based algorithms and an Horvitz-Thompson estimator to assess the heterogeneity of treatment and spillover effects with respect to individual, neighborhood and network characteristics in the context of clustered network interference. We illustrate how the proposed binary tree methodology performs in a Monte Carlo simulation study. Additionally, we provide an application on a randomized experiment aimed at assessing the heterogeneous effects of information sessions on the uptake of a new weather insurance policy in rural China.

## Getting Started

Installing the latest developing version: 

```r
library(devtools)
install_github("fbargaglistoffi/Network-Causal-Tree", ref="master")
```

Import:

```r
library("Network-Causal-Tree")
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
