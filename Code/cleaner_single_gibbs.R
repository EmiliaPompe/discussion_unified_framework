## This script aims at simply running the Gibbs sampler on a small data set
rm(list = ls())
set.seed(2020)
library(Rcpp)
library(tidyverse)
theme_set(theme_minimal())
source("getData.R")
sourceCpp("init_clustering.cpp")
sourceCpp("compute_loglikelihood_clusters.cpp")
sourceCpp("update_eta.cpp")
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))
source("relabel.R")
##


## The data is the matrix V
## important for indexing: entries of V should start at zero
V <- V - 1 
## let's work with a small data set for the moment
V <- V[1:10,1:4]
## numbers of possibilities for each column
## this is denoted by M_ell in the overleaf document
dimV <- as.integer(dimV[1:4])
## the following vector is...
cumdime <- 1 +  c(0,  cumsum(dimV)) 
## define dimension of V
n <- dim(V)[1]
p <- dim(V)[2]
## 
## prior parameter on N 
g <- 1.02
##
## The state of the Markov chain  is (eta, N, b', b0, theta)
## initial eta 
eta <- sample(n, size = n, replace = T) - 1
partition <- init_clustering(eta)
relabel_result <- relabel2(eta, partition)
eta <- relabel_result$eta
partition$clsize <- partition$clsize[relabel_result$old_to_new] 
partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]

# clusteval::cluster_similarity(eta, relabel_result$eta)
# identical(init_clustering(relabel_result$eta), relabel_result$clustering)
## initial N
N <- 50
m0 <- logit(0.01)
s0 <- sqrt(0.1)
## initial b0
beta0 <- rnorm(p, m0, s0)
s  <- sqrt(0.5)
## initial beta
beta <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p) beta[,field] <- rnorm(n, beta0[field], s)
alpha <- expit(beta)
## initial frequencies of categories, initial values
theta <- as.double(unlist(ALPHA[1:p])) 
## compute log-likelihood associated with each cluster
partition_ll <- compute_loglikelihood_clusters(partition, theta, V, dimV, alpha)

## next iterate updates in the Gibbs sampler
## update of lambda
update_eta_result <- update_eta(eta, partition, partition_ll, theta, V, dimV, alpha, N)
eta <- update_eta_result$eta
partition_ll <- update_eta_result$clusterloglikelihoods
partition <- update_eta_result$clustering
## relabel
relabel_result <- relabel2(eta, partition)
eta <- relabel_result$eta
partition$clsize <- partition$clsize[relabel_result$old_to_new] 
partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
partition_ll <- partition_ll[relabel_result$old_to_new,]
## update of N
# random walk proposal
delta <- 10
Nproposal <- sample(x = (N-delta):(N+delta), size = 1)


## update of beta 

alpha <- expit(beta)
## update of beta0 



