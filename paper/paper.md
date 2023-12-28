---
title: 'NetworkCausalTree: An R package for heterogeneous treatment effects under network interference'
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
    affiliation: "3"
affiliations:
 - name: Sant'Anna School for Advanced Studies
   index: 1
 - name: Department of Biostatistics, Harvard School of Public Health
   index: 2
 - name: Causal Learning and Artificial Intelligence lab, Institute of Science and Technology Austria
   index: 3
 - name: Department of Biostatistics, Yale School of Public Health
   index: 4
date: 01 January 2024
bibliography: paper.bib
---

# Summary

@cox1958planning asserts that interference occurs when the assignment of treatment to one unit influences the outcomes of other units. In the realm of policy interventions, interference can manifest through various interactions, encompassing social, physical, or virtual connections. The conventional Rubin Causal Model, employed in causal inference studies [@rubin1986comment], excludes interference. However, when interference is present, it must be taken into account to prevent a biased estimation of the treatment effect [@forastiere2016identification] and to be able to estimate spillover effects, that is, the effects on a unit's outcome of the treatment received by other units. Consequently, recent research has devised new methodologies to tackle interference on networks [@forastiere2016identification; @aronow2017estimating; @athey2018exact; @leung2020treatment].
Concurrently, alongside the interference research domain, scholars have crafted machine learning algorithms to appraise treatment effects' heterogeneity with respect to individual characteristics [@athey2016recursive]. These algorithms identify sub-populations by iteratively segregating groups whose estimated average treatment effect deviates the most. Understanding these heterogeneities enables researchers to understand which sub-populations are particularly responsive to the intervention and, consequently, to improve the design of future targeting strategies.

To integrate the aforementioned two topics in the field of causal inference, [@bargagli2020heterogeneous] introduces a novel machine learning algorithm, named Network Causal Tree (NCT), that explores the heterogeneity of treatment and spillover effects concerning individual, neighborhood, and network characteristics within randomized settings. Evaluating the heterogeneity of treatment and spillover effect is crucial to identify both those sub-populations who are responsive to the intervention they actually received and those who are highly susceptible to the treatment received by other units. Heterogeneity exploration is performed via an interpretable decision tree algorithm. Decision trees are _interpretable_ in the sense proposed by [@murdoch2019definitions] as they provide _readily available information on the relationships that they have learned from the data_. NCT is designed to operate within clustered network interference (CNI), where agents are partioned into distinct clusters, and spillover mechanisms exclusively take place within clusters based on the links of a cluster-specific network. The estimation of conditional effects is carried out using an extended version of the Horvitz-Thompson estimator [@aronow2017estimating], tailored to accommodate clustered network interference. `NetworkCausalTree` is an R Package providing a flexible implementation of the Network Causal Tree algorithm. 

# Statement of need

The existing R algorithms are designed either to detect and estimate Heterogeneous Treatment Effects (HTEs) or Heterogeneous Spillover Effects (HSEs) in the presence of clustered network interference.

On the one hand, there are a plethora of methods to discover and estimate HTEs employing decision trees in the absence of intereference. Some of them have also been implemented in R. For instance: i) the `aVirtualTwins` package, which implements the Virtual Twins algorithm introduced by [@foster2011subgroup]; ii) the `causalTree` package, which implements the Causal Honest Tree algorithm introduced by Athey and Imbens [@athey2016recursive]; iii) the `policytree` which implements the Policy Tree methodology as proposed by [@athey2021policy]; iv) the `CRE` package [@Cadei2023], which implements the Causal Rule Ensemble algorithm introduced by [@bargagli2020causal]. 

On the other hand, there are two main packages to estimate the HTEs and HTEs in the presence of interference: i) the `inferference` package implements the Inverse Probability Weighting (IPW) estimators proposed by [@tchetgen2012causal] to evaluate the causal effect of an intervention in the presence of partial interference; ii) the `interference` package estimates causal effects under network interference [@aronow2017estimating]. 

To account for the current lack of existing methods for a data-driven discovery of the heterogeneity of causal effects in the presence of network interference, we propose `NetworkCausalTree`, an R package providing a flexible implementation of the Network Causal Tree algorithm. `NetworkCausalTree` provides i) an interpretable representation of the heterogeneity of treatment and spillover effects in the presence of CNI, within randomized settings; ii) the possibility for the researcher to simultaneously model the heterogeneity of both treatment and spillover effects. To the best of our knowledge, there aren't existing packages to explore, in a data-driven manner, the heterogeneity of treatment and spillover effects in the presence of network interference. 

# Algotithm 

The NCT algorithm is designed to detect and estimate heterogeneous treatment and spillover effects in randomized settings under CNI. Under CNI, spillover is confined to units within the same cluster, and an individual's outcome is influenced by the treatment status of directly connected units in the cluster-specific network. Potential outcomes $ Y_{ik}(w,g)$ - where $k\in \mathcal{K}=[1, \ldots, K]$ is the cluster indicator, and $i=1,\ldots, n_k$ is the unit indicator in each cluster $k$ - are indexed with respect to the individual intervention $W_{ik}=w$ and to the neighborhood treatment $G_{ik}=g$, which is defined as a binary variable that equals 1 if the unit $ik$ has at least one treated neighbor, 0 otherwise. For each observation $ik$ in the sample $\mathcal{V}$ we observe a vector of individual or network characteristics $\mathbf{X_{ik}}$.  The adjacency matrix describing interactions within a given cluster $k$ is denoted with $\mathbf{A_{k}}$. Estimands of interest ($\tau_{(w,g;w',g')}$) encompass two treatment and two spillover effects, comparing average potential outcomes under varying individual or neighborhood treatment statuses, while keeping the other constant. 

$$
\tau_{(w,g;w',g')}(l(\mathbf{x})) =  \frac{1}{N(l(\mathbf{x}))}\sum_{k=1}^K\sum_{i=1}^{n_k}Y_{ik}(w,g)\mathds{1}(\mathbf{X}_{ik}\in l(\mathbf{x)) -  \frac{1}{N(l(\mathbf{x}))}\sum_{k=1}^K\sum_{i=1}^{n_k}Y_{ik}(w',g')\mathds{1}(\mathbf{X}_{ik} \in l(\mathbf{x})). 
$$

These effects will be estimated within each leaf through an appropriate Horvits-Thompson estimator [@aronow2017estimating], which employs both the marginal probability of an individual $ik$ to be exposed to the treatment level $(w,g)$ $\pi_{ik}(w,g)$ and the joint probablity of two individuals $ik$ and $jk'$ to be simultaneously exposed to the treatment levels $(w,g)$ and $(w',g')$, respectively, $\pi_{ikjk'}(w,g;w',g')$. 

`NetworkCausalTree` relies on a flexible in-sample splitting criterion function that can potentially accommodate more than one estimator simultaneously. `NetworkCausalTree` procedure is divided into two steps: discovery and estimation. Each cluster is randomly assigned to either the discovery or the estimation set (using a procedure that is commonly referred to as _honest splitting_). In the discovery step, `NetworkCausalTree` builds a tree according to the in-sample splitting criterion and stops when either the tree has reached its maximum depth or any additional split would generate leaves, which are not sufficiently representative of the four exposure conditions. In the estimation step, `NetworkCausalTree` estimates the CACEs and their standard errors in each leaf. The steps of the `NetworkCausalTree` algorithm are illustrated in the Algorithm 1 table.

![Illustration of the Network Causal Tree Algorithm.](images/NCT_algorithm.png)

# Usage

`NetworkCausalTree` is available on [GitHub](https://github.com/fbargaglistoffi/NetworkCausalTree) and can be installed and loaded into the R session
using:


```r
library(devtools)
install_github("fbargaglistoffi/NetworkCausalTree", ref="master")
```

`data_generator()` is a flexible synthetic dataset generator, which can be used for simulations before applying `NetworkCausalTree` to real-world  data sets. It returns a CNI environment, where cluster-specific networks are generated either using the Erdos-Renyi model [@erdos1959random], the Barabasi-Albert [@barabasi1999emergence], or the Exponential Random Graph Model[@lusher2013exponential], that can potentially account for homophily.

```r
dataset <- data_generator(N = 4000, 
                          K = 4,
                          m = 80, 
                          p = rep(0.2,4000), 
                          het = TRUE, 
                          h = 2, 
                          method_networks = "er", 
                          param_er = 0.1)
```

We propose here two examples of how to run the Network Causal Tree algorithm by the `NetworkCausalTree` package.


**Example 1.** Running Network Causal Tree while relying on a singular in-sample splitting criterion function, that assesses the heterogeneity of the main treatment effect ($\tau_{(1,0;0,0)}$)  only 


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

**Example 2.** Running Network Causal Tree while relying on a composite in-sample splitting criterion function, which accounts for all the four estimands simultaneously 

```r
result <- NetworkCausalTree(X = dataset[["X"]],
                            Y = dataset[["Y"]],
                            W = dataset[["W"]], 
                            A = dataset[["A"]],
                            M = dataset[["M"]],
                            p = dataset[["p"]], 
                            effect_weights = c(0.25,0.25,0.25,0.25),
                            ratio_dis = 0.5,
                            depth = 3,
                            minsize = 5, 
                            method = "composite",
                            output = "estimation")
```


The results are included in a `data.frame` object which provides information about i) the partitions identified by the tree, ii) the number of units included in all partitions, iii) the estimated CACEs, with their corresponding standard errors, in all leaves. The function `print()` displays these results and the `plot_NCT()` function visualizes results. For instance, the plot related to the `NetworkCausalTree` resulted in the **Example 1** is obtained using: 

```r
title <- expression(paste("CAUSAL TREE TARGETED TO ",tau,"(1,0;0,0)"),sep="")
cov_names <- colnames(dataset[["X"]])

plot_NCT(NCT = result, 
         cov_names = cov_names,
         title = title)
```

Figure 1 reports the visualization of the results for Example 1, which discovers correctly the underlying CATE decomposition.

![Network Causal Tree targeted to the main treatment effect](images/result.jpeg)

Online documentation for the package can be found at [https://github.com/fbargaglistoffi/NetworkCausalTree](https://github.com/fbargaglistoffi/NetworkCausalTree).

# Acknowledgements

This work was supported by the Alfred P. Sloan Foundation Grant G-2020-13946 for the development of _Causal Inference with Complex Treatment Regimes: Design, Identification, Estimation, and Heterogeneity_ and the Harvard Data Science Initiative Postdoctoral Research Fund Award.

# References
