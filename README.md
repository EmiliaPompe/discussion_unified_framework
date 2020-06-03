## *[A Unified Framework for De-Duplication and Population Size Estimation](https://projecteuclid.org/euclid.ba/1551949260)* 

### A discussion by [Nianqiao Ju](https://phylliswithdata.com/about/), Niloy Biswas, [Pierre E. Jacob](https://sites.google.com/site/pierrejacob/), [Gonzalo Mena](http://gomena.github.io), John O’Leary and [Emilia Pompe](https://www.stats.ox.ac.uk/~pompe/).


This repository contains documentation and code for the discussion paper [Unified Framework for De-Duplication and Population Size Estimation](https://projecteuclid.org/euclid.ba/1551949260)[1], authored by **Andrea Tancredi**, **Rebecca Steorts**  and **Brunero Liseo**,  to appear in **Bayesian Analysis**. It is organized as follows.

### Documentation
We provide two write-ups in the [`document`](https://github.com/EmiliaPompe/discussion_unified_framework/tree/master/document) directory. First, the [comment on the discussion paper](https://github.com/EmiliaPompe/discussion_unified_framework/blob/master/document/badiscussion.pdf), 
in which we state the possibility of using recently developed methodology [2,3], based on coupling of two chains, applied to the Gibbs sampler proposed in Section 5 of [[1]](https://projecteuclid.org/euclid.ba/1551949260). Then, in a [second document](https://github.com/EmiliaPompe/discussion_unified_framework/blob/master/document/documentation.pdf), we give more details on the Gibbs sampler, a coupling of it, and we present some empirical results. These indicate that coupling techniques are applicable, as a proof of concept, and can provide useful convergence assesments following [2].

### Code
Code reproducing our experiments and presented in the format of a `R` package is contained in the [`couplingdeduplication`](https://github.com/EmiliaPompe/discussion_unified_framework/tree/master/couplingdeduplication) directory. These pertain the synthetic dataset `RLdata500` from the `R` package `RecordLinkage` analysed in Section 6 of [[1]](https://projecteuclid.org/euclid.ba/1551949260).

Some explanations:

- inst/run_authorscode.R: runs the authors', either keeping both theta and beta fixed, or just beta, and saves the output 

- inst/run_gibbs.R: runs this implementation, either keeping both theta and beta fixed, or just beta, and saves the output 

- inst/compare_gibbs.R: creates plots to check the similarity of results between this implementation and the authors', provided
one has executed 'run_gibbs.R' and 'run_authorscode.R' already

- inst/run_coupledgibbs.R: runs pairs of Markov chains until they meet, in order to obtain meeting times and TV upper bounds

- inst/run_coupledgibbs_norelabel.R: runs pairs of Markov chains until they meet, without the relabeling step...
to see that indeed the relabeling step is useful.


### References
[1]Andrea Tancredi, Rebecca Steorts, and Brunero Liseo. [A Unified Framework for De-Duplication and Population Size Estimation Bayesian Analylsis, advance publication](https://projecteuclid.org/euclid.ba/1551949260), 2020.

[2] Niloy Biswas, Pierre E Jacob, and Paul Vanetti. [Estimating convergence of Markov chains with L-lag couplings. In Advances in Neural Information Processing Systems](https://papers.nips.cc/paper/8958-estimating-convergence-of-markov-chains-with-l-lag-couplings), pages 7389–7399, 2019. 

[3] Pierre E Jacob, John O’Leary, and Yves F Atchad ́e. [Unbiased Markov chain Monte Carlo with couplings. Journal of the Royal Statistical Society: Series B (Statistical Methodology) (with discussion) (to appear)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12336), 2020. 








