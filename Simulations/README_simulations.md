# Monte Carlo Simulations

The codes <a href="https://github.com/fbargaglistoffi/Network-Causal-Tree/tree/master/Simulations">in this folder</a> implement a number of Monte Carlo simulations to investigate the fit of the _Network Causal Tree_ algorithm.

The general data generating process (dgp) is described in detail <a href="https://github.com/fbargaglistoffi/Network-Causal-Tree/tree/master/Simulations">in the paper</a>. In a nutshell, for each simulation we build a number of clusters each of size 100, generating an Erdős–Rényi network within each cluster, and 10 covariates. </br>
Out of these generated clusters, a half are used for the _discovery sample_ and the other half for the _estimation sample_.
We introduce two and four _true causal rules_, each of them with depth two (e.g.,  _x1==0 & x2==0_, namely _x1_ and _x2_ are the two _heterogeneity driving variables_ [HDVs]).

The results are reported in the paper, and were obtained by aggregating the results over 500 different datasets for each effect size.

## R code

The code files used to implement the simulations are the following:
* <tt>`simulations_main_spillover`</tt>:  two causal rules for the main and spillover effect with two different HDVs (_x1, x2_), uncorrelated regressors, 1000/2000/3000 data points and 10/20/30 clusters;
* <tt>`simulations_main_spillover_correlated`</tt>: two causal rules for the main and spillover effect with two different HDVs (_x1, x2_), correlated regressors (<tt>`rho`</tt> = 0.25/0.50), 3000 data points and 30 clusters;
* <tt>`simulations_main_homophily`</tt>:  two causal rules for the main and spillover effect with two different HDVs (_x1, x2_), homophily network, uncorrelated regressors, 3000 data points and 30 clusters;
* <tt>`simulations_different_main_spillover`</tt>: four causal rules for the main and spillover effect with three different HDVs (_x1, x2, x3_), uncorrelated regressors, 3000 data points and 30 clusters.

### Parameters used for the simulations

* <tt>`N`</tt>: number of data points (1000/2000/3000);
* <tt>`M`</tt>: number of clusters (10/20/30);
* <tt>`n_cov`</tt>: number of covariates (10);
* <tt>`rho`</tt>: correlation within the covariates (0/0.25/0.50);
* <tt>`prob`</tt>: treatment assignment probability (0.50);
* <tt>`crls`</tt>: number of causal rules (2/4);
* <tt>`seq`</tt>: effect size magnitude (varying between 0.1 and 10.1);
* <tt>`nsim`</tt>: number of datasets created (500);
* <tt>`gsize`</tt>: size of each cluster (100).

