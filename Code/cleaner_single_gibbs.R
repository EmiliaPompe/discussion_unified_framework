## This script aims at simply running the Gibbs sampler on a small data set
## this is work in progress
rm(list = ls())
set.seed(1)
## load packages
library(Rcpp)
library(tidyverse)
## set graphical themes
theme_set(theme_minimal())
## load data
source("getData.R")
## source functions
## to relabel 
source("relabel.R")
## to define clustering/partition given vector of labels 'eta'
sourceCpp("init_clustering.cpp")
## to compute log-likelihood of the different clusters in the partition
sourceCpp("compute_loglikelihood_clusters.cpp")
## to perform full sweep of eta updates 
sourceCpp("update_eta.cpp")
## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))
##
## The data is the matrix V
## important for indexing: entries of V should start at zero, following C convention
V <- V - 1 
## let's work with a small data set for the moment
V <- V[1:50,1:2]
## numbers of possibilities for each column; denoted by M_ell in the overleaf document
dimV <- as.integer(dimV[1:ncol(V)])
## define dimensions of V
n <- dim(V)[1]
p <- dim(V)[2]
## prior parameter for N 
g <- 1.02
## number of MCMC iterations
nmcmc <- 1e4
## whether to print some things during the run, or not
verbose <- TRUE
## record history of certain components
N_history <- rep(NA, nmcmc)
theta1_history <- matrix(NA, nrow = nmcmc, ncol = dimV[1])
##
## The state of the Markov chain  is (eta, N, b', b0, theta)
## First we initialize the chains
## initial N
## initial eta 
eta <- sample(x = 1:n, size = n, replace = T) - 1
N <- 2 * length(unique(eta))
partition <- init_clustering(eta)
relabel_result <- relabel2(eta, partition)
eta <- relabel_result$eta
partition$clsize <- partition$clsize[relabel_result$old_to_new] 
partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
# clusteval::cluster_similarity(eta, relabel_result$eta)
# identical(init_clustering(relabel_result$eta), relabel_result$clustering)
m0 <- logit(0.01)
s0 <- sqrt(0.1)
## initial b0
beta0 <- rnorm(p, m0, s0)
s  <- sqrt(0.5)
## initial beta
beta <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p) beta[,field] <- rnorm(n, beta0[field], s)
alpha <- expit(beta)
alpha[which(partition$clsize==0),] <- -Inf
## initial frequencies of categories, initial values
theta <- ALPHA[1:p]
## compute log-likelihood associated with each cluster
partition_ll <- compute_loglikelihood_clusters(partition, theta, V, dimV, alpha)
## initialize quantities to monitor
Naccept <- 0
thetaaccept <- rep(0, p)
## next iterate updates in the Gibbs sampler
for (imcmc in 1:nmcmc){
  if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, "\n")
  ## update of eta
  # print(N)
  # print(length(unique(eta)))
  update_eta_result <- update_eta(eta, partition, partition_ll, theta, V, dimV, alpha, N, beta0, s)
  eta <- update_eta_result$eta
  partition_ll <- update_eta_result$clusterloglikelihoods
  partition <- update_eta_result$clustering
  alpha <- update_eta_result$alpha
  ## relabel 
  relabel_result <- relabel2(eta, partition)
  eta <- relabel_result$eta
  partition$clsize <- partition$clsize[relabel_result$old_to_new]
  partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
  partition_ll <- partition_ll[relabel_result$old_to_new,]
  alpha <- alpha[relabel_result$old_to_new,]
  ## update of N
  ## random walk proposal
  N_rw_stepsize <- 20
  # if (verbose) print("update N")
  Nproposal <- sample(x = (N-N_rw_stepsize):(N+N_rw_stepsize), size = 1)
  if (Nproposal >= partition$ksize){
    target_proposal <- lfactorial(Nproposal) - lfactorial(Nproposal - partition$ksize) - (n+g) * log(Nproposal)
    target_current <- lfactorial(N) - lfactorial(N - partition$ksize) - (n+g) * log(N)
    logu <- log(runif(1))
    if (logu < (target_proposal - target_current)){
      N <- Nproposal
      Naccept <- Naccept + 1
    }
  }
  N_history[imcmc] <- N
  # if (verbose) print("update N done")
  ## To add: update of beta 
  ##
  ## alpha <- expit(beta)
  ## partition_ll <- compute_loglikelihood_clusters(partition, theta, V, dimV, alpha)
  ## To add: update of beta0 
  
  ## update of theta
  ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
  ## but for numerical stability we remove the part close to the boundary 
  for (field in 1:p){
    ## theta associated with field
    theta_current <- theta[[field]]
    ## propose modification using Dirichlet proposal
    theta_scale <- 450
    theta_proposal <- rgamma(n = dimV[field], shape = theta_current * theta_scale, rate = 1)
    theta_proposal <- theta_proposal / sum(theta_proposal)
    ## compute Dirichlet log density
    dirichlet1 <- function(x, alpha) {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s <- (alpha - 1) * log(x)
      s <- ifelse(alpha == 1 & x == 0, -Inf, s)
      return(sum(s) - logD)
    }
    ## compute MH ratio
    if (any((theta_proposal < 1e-10) | (theta_proposal > (1-1e-10)))){
      ## reject move as too close to boundary
      # print("reject because too close to boundary")
      # print(theta_current)
      # print(theta_proposal)
    } else {
      proposal_logratio <- dirichlet1(theta_current, theta_proposal * theta_scale) - dirichlet1(theta_proposal, theta_current * theta_scale)
      # print(proposal_logratio)
      partition_ll_field_proposal <- compute_loglikelihood_onefield(field-1, partition, theta_proposal, V, dimV, alpha)
      mh_ratio <- sum(partition_ll_field_proposal[partition$clsize>0]) - sum((partition_ll[,field])[partition$clsize>0]) + proposal_logratio
      ## accept proposal or not
      if (log(runif(1)) < mh_ratio){
        theta[[field]] <- theta_proposal
        partition_ll[,field] <- partition_ll_field_proposal
        thetaaccept[field] <- thetaaccept[field] + 1
      }
    }
  }
  theta1_history[imcmc,] <- theta[[1]]
  #
}
 
cat(Naccept/nmcmc, "\n")
cat(thetaaccept/nmcmc, "\n")

matplot(N_history, type = 'l')

matplot(theta1_history[,1:2], type = 'l')
abline(h = ALPHA[[1]][1:2])

# pairs(theta1_history[200:nmcmc,1:10])
# hist(N_history[(nmcmc/10):nmcmc])

hist(theta1_history[(nmcmc/10):nmcmc,1])
abline(v = ALPHA[[1]][1])

hist(theta1_history[(nmcmc/10):nmcmc,2])
abline(v = ALPHA[[1]][2])

hist(theta1_history[(nmcmc/10):nmcmc,3])
abline(v = ALPHA[[1]][3])

 