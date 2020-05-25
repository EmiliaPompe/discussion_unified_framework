## *[A Unified Framework for De-Duplication and Population Size Estimation](https://projecteuclid.org/euclid.ba/1551949260)*: a discussion.
 
 #### by [Nianqiao Ju](https://phylliswithdata.com/about/), Niloy Biswas, [Pierre E. Jacob](https://sites.google.com/site/pierrejacob/), [Gonzalo Mena](http://gomena.github.io), John O’Leary and [Emilia Pompe](https://www.stats.ox.ac.uk/~pompe/).


This repository contains documentation and code for the discussion paper [Unified Framework for De-Duplication and Population Size Estimation](https://projecteuclid.org/euclid.ba/1551949260), authored by **Andrea Tancredi**, **Rebecca Steorts**  and **Brunero Liseo**,  to appear at **Bayesian Analysis**. It is organized at follows

### Documentation
We provide two write-ups [/document](https://github.com/EmiliaPompe/discussion_unified_framework/tree/master/document). First, in the [official discussion paper](https://github.com/EmiliaPompe/discussion_unified_framework/blob/master/document/badiscussion.pdf ), we state the possibility of using the recently developed coupling [1,2] methodology for chains produced with the Gibbs sampler introduced in section 5 of the [original paper](https://projecteuclid.org/euclid.ba/1551949260). Then, on a [second document](https://github.com/EmiliaPompe/discussion_unified_framework/blob/master/document/documentation.pdf), we give a detailed exposition of such sampler on a simple setup, and show empirical results through numerical experiments. These results indicate that the presented coupling leads to similar results as the one presented in the original paper (Figure 1), where convergence is assessed according to the novel diagnostics methodology described in [1] (Figure 2).

### Code
Code reproducing our experiments and presented as a R package is contained in [/couplingdeduplication](https://github.com/EmiliaPompe/discussion_unified_framework/tree/master/couplingdeduplication).

### References

[1] Niloy Biswas, Pierre E Jacob, and Paul Vanetti. Estimating convergence of Markov chains with L-lag couplings. In Advances in Neural Information Pro- cessing Systems, pages 7389–7399, 2019. 

[2] Pierre E Jacob, John O’Leary, and Yves F Atchad ́e. Unbiased Markov chain Monte Carlo with couplings. Journal of the Royal Statistical Society: Series B (Statistical Methodology) (with discussion) (to appear), 2020. 








