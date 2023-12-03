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

Summary of the paper.

# Statement of need

Why `NetworkCausalTree` package is necessary.

# Usage

[@bargagli2020heterogeneous]


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
---
bibliography: paper.bib
nocite: "@*"
---

# Acknowledgements

We acknowledge ...

# References
