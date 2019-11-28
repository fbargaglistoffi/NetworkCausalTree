# Monte Carlo Simulations

This code implements a number of Monte Carlo simulations to investigate the fit of the _Network Causal Tree_ algorithm.
The general data generating process (dgp) is the following: for each simulation we built 15 clusters (of size 100 and 200, respectively), generating an Erdős–Rényi network within each cluster, and 10 covariates.
Out of these 15 clusters, 5 were used for the _discovery sample_, 5 for the _estimation sample_ and 5 for the _test sample_.
We introduced two _true causal rules_, each of them with depth two (e.g.,  x1==0 & x2==0, namely _x1_ and _x2_ are the two _heterogeneity driving variables_ [HDVs]).
The results are reported in a number of _excel tables_, and were obtained by aggregating the results, for each effect size, on 100 different datasets.

## Parameters used for the simulations

* <tt>`N`</tt>: number of data points (1500/3000)
* <tt>`n_cov`</tt>: number of variables (10)
* <tt>`rho`</tt>: correlation within the covariates (0/0.25)
* <tt>`prob`</tt>: treatment assignment probability (0.4)
* <tt>`crls`</tt>: number of causal rules (2/3)
* <tt>`seq`</tt>: effect size magnitude (varying between 0.1 and 10.1)
* <tt>`nsim`</tt>: number of datasets created (100)
* <tt>`gsize`</tt>: size of each cluster (100/200)

## Results

The simualations' results were obtained for both sample sizes 1500 and 3000.
Below a brief guide to the _excel files_:
* <tt>`one_main_effect`</tt>: two causal rule for just one main effect (_tau(10,00)_) with two different HDVs (_x1, x2_)
* <tt>`two_main_effects`</tt>: two causal rules for the two main effects (_tau(10,00)_, _tau(11,01)_) with two different HDVs (_x1, x2_)
* <tt>`two_main_effects_correlated`</tt>: two causal rules for the two main effects (_tau(10,00)_, _tau(11,01)_) with two different HDVs (_x1, x2_) and correlated regressors (<tt>`rho`</tt> = 0.25)
* <tt>`three_main_effects_overlap`</tt>: three causal rules for the two main effects (_tau(10,00)_, _tau(11,01)_), with three different HDVs (_x1, x2, x3_)
* <tt>`three_main_effects_overlap`</tt>: three causal rules for the two main effects (_tau(10,00)_, _tau(11,01)_), with three different HDVs (_x1, x2, x3_), and different effect sizes (the effects for _tau(11,01)_ are 50% bigger than the effects for _tau(10,00)_)
* <tt>`two_spillover_effects`</tt>: two causal rules for the two spillover effects (_eta(11,10)_, _eta(01,00)_) with two different HDVs (_x1, x2_)

The following results were reported:
* <tt>`correct_rules`</tt>: number of rules correctly detected by the algorithm;
* <tt>`mse_tau_est_1000`</tt>: Monte Carlo mean-squared-error for _tau(10,00)_ in the estimation sample
* <tt>`bias_tau_est_1000`</tt>: Monte Carlo bias for _tau(10,00)_ in the estimation sample
* <tt>`mse_tau_test_1000`</tt>: Monte Carlo mean-squared-error for _tau(10,00)_ in the test sample
* <tt>`bias_tau_test_1000`</tt>: Monte Carlo bias for _tau(10,00)_ in the test sample
* <tt>`mse_tau_est_1101`</tt>: Monte Carlo mean-squared-error for _tau(11,01)_ in the estimation sample
* <tt>`bias_tau_est_1101`</tt>: Monte Carlo bias for _tau(11,01)_ in the estimation sample
* <tt>`mse_tau_test_1101`</tt>: Monte Carlo mean-squared-error for _tau(11,01)_ in the test sample
* <tt>`bias_tau_test_1101`</tt>: Monte Carlo bias for _tau(11,01)_ in the test sample
* <tt>`coverage_est_1000`</tt>: Monte Carlo coverage for _tau(10,00)_ in the estimation sample
* <tt>`coverage_test_1000`</tt>: Monte Carlo coverage for _tau(10,00)_ in the test sample
* <tt>`coverage_est_1101`</tt>: Monte Carlo coverage for _tau(11,01)_ in the estimation sample
* <tt>`coverage_test_1101`</tt>: Monte Carlo coverage for _tau(11,01)_ in the test sample
* <tt>`vi_x1`</tt>: variable importance for _x1_ 
* <tt>`vi_x2`</tt>: variable importance for _x2_ 
* <tt>`vi_X`</tt>: variable importance for all the variables with the exclusion of _x1_ and _x2_ 

N.B.: (i) the same results are reported for spillover effect but instead of _tau(10,00)_, _tau(11,01)_ we have _eta(11,10)_, _eta(01,00)_; (ii) the variable importance refers to the average increase in splitting criterion due to split on that variable (the larger the higher the importance of the variable); (iii) in case of <tt>`three_main_effects_overlap`</tt> the variable importance is reported also for _x3_.





