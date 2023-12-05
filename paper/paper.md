---
title: 'NetworkCausalTree: an R package for ...'
tags:
  - R
  - TODO
authors:
  - name: Costanza Tortù
    orcid: 0000-0003-2561-9726
    equal-contrib: true
    affiliation: "1"
  - name: Falco J. Bargagli Stoffi
    orcid: 0000-0002-6131-8165
    equal-contrib: true
    affiliation: "2"
  - name: Riccardo Cadei
    orcid: 0000-0003-2416-8943
    affiliation: "2"
  - name: Laura Forastiere
    orcid: 0000-0003-3721-9826
    affiliation: "2"
affiliations:
 - name: Sant'Anna School for Advanced Studies
   index: 1
 - name: Department of Biostatistics, Harvard School of Public Health
   index: 2
date: 24 March 2023
bibliography: paper.bib
---

# Summary

[@cox1958planning] asserts that interference occurs when the assignment of treatment to one unit influences the outcomes of other units. In the realm of policy interventions, interference can manifest through various interactions, encompassing social, physical, or virtual connections. The conventional Rubin Causal Model, employed in causal inference studies [@rubin1986comment], excludes interference. However, when interference plays a role, it introduces bias into estimates [@forastiere2016identification]. Additionally, spillover effects enable researchers to gauge the overall impact of an intervention and to enhance the efficiency of treatment assignment mechanisms. Consequently, recent research has devised inventive methodologies to tackle interference [@sobel2006randomized,@rosenbaum2007interference].
Concurrently, alongside the interference research domain, scholars have crafted machine learning algorithms to appraise treatment effects' heterogeneity with respect of individual characteristics [@athey2016recursive]. The rationale behind these algorithms lies in partitioning sub-populations by iteratively segregating groups whose estimated average treatment effect exhibits the most deviation.

To integrate the aforementioned two topics in the field of causal inference through [@bargagli2020heterogeneous] introduces a novel machine learning algorithm, named Network Causal Tree (NCT), that explores the heterogeneity of treatment and spillover effects concerning individual, neighborhood, and network characteristics within randomized settings. NCT is designed to operate amidst clustered network interference, where agents are categorized into distinct clusters, and spillover mechanisms exclusively take place within clusters based on the links of a cluster-specific network. The estimation of conditional effects is carried out using an extended version of the Horvitz-Thompson estimator [@aronow2017estimating], tailored to accommodate clustered network interference. `NetworkCausalTree` is an R Package providing a flexible implementation of the Network Causal Tree algorithm.


# Statement of need

Why `NetworkCausalTree` package is necessary.

# Usage




```r
library(devtools)
install_github("fbargaglistoffi/NetworkCausalTree", ref="master")
```

```r
dataset <- data_generator(N = 4000, 
                          K = 4,
                          m = 80, 
                          p = rep(0.2,2000), 
                          het = TRUE, 
                          h = 2, 
                          method_networks = "er", 
                          param_er = 0.1)
```

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
```

```r
title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
cov_names <- colnames(dataset[["X"]])

plot_NCT(NCT = result, 
         cov_names = cov_names,
         title = title)
```


# Acknowledgements

We acknowledge ...

# References
