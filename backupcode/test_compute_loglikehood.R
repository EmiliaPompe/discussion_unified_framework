## This script aims at simply running the Gibbs sampler on a small data set
## this is work in progress
rm(list = ls())
set.seed(1)
## load packages
library(Rcpp)
library(tidyverse)
## set graphical themes
theme_set(theme_minimal())
## load data 'V' and defines n and p, the dimensions of V
library(RecordLinkage)
data(RLdata500)
myRLDATA=cbind(id=identity.RLdata500,RLdata500)
M1=matrix(unlist(strsplit(soundex(myRLDATA$fname_c1),split="")),ncol=4,byrow=TRUE)
M2=matrix(unlist(strsplit(soundex(myRLDATA$lname_c1),split="")),ncol=4,byrow=TRUE)
M3=matrix(unlist(strsplit(as.character(myRLDATA[,6]),split="")),ncol=4,byrow=TRUE)
M4=matrix( character(nrow(myRLDATA)*2),ncol=2)
for (i in 1:2) M4[,i]=as.character(myRLDATA[,6+i])

V=cbind(M1,M2,M3,M4)
V=as.data.frame(V)
p=ncol(V);

fieldfrequencies=list()
for (l in 1:p) fieldfrequencies[[l]]=table(V[,l])/nrow(V)
Mvec <- c()
for (l in 1:p) {
  Mvec[l]=length(levels(V[,l]))
  V[,l]=as.numeric(V[,l])
}
V=as.matrix(V)
# V is a matrix of n * p
# entries are values in 1:M_l, where M_l is the number of possible categories in field l
# and field l in in 1:p
# Vector 'Mvec' has entries 'M_l' for l in 1:p

## let's subset the data 
V <- V[1:50,1:2]
## define dimensions of V
n <- dim(V)[1]
p <- dim(V)[2]

## get some 'ground truth'
# MATCH=outer(myRLDATA$id, myRLDATA$id, FUN="==")
# MATCH=MATCH[upper.tri(MATCH)]

## source functions
## to relabel 
source("relabel.R")
## to define clustering/partition given vector of labels 'eta'
sourceCpp("init_clustering.cpp")
## to compute log-likelihood of the different clusters in the partition
sourceCpp("compute_loglikelihood_all_clusters_all_fields.cpp")
sourceCpp("compute_loglikelihood_all_clusters_one_field.cpp")
sourceCpp("compute_loglikelihood_one_cluster_one_field.cpp")
## to perform full sweep of eta updates 
sourceCpp("update_eta.cpp")
## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))


## hyper parameter specification
## prior parameter for N 
g <- 1.02
## prior parameter for beta0
m0 <- logit(0.01)
s0 <- sqrt(0.1)
## prior parameter for beta given beta0
s  <- sqrt(0.5)

## number of MCMC iterations
nmcmc <- 1e4
## there should be update frequencies ... 

## whether to print some things during the run, or not
verbose <- TRUE

## record history of certain components
N_history <- rep(NA, nmcmc)
theta1_history <- matrix(NA, nrow = nmcmc, ncol = Mvec[1])

## The state of the Markov chain  is (eta, N, b', b0, theta)

## Initialization of the chains
## initial eta 
eta <- sample(x = 1:n, size = n, replace = T)
## initial N
N <- max(n, 2 * length(unique(eta)))
## 
partition <- init_clustering_cpp(eta-1)
relabel_result <- relabel2(eta, partition)
# clusteval::cluster_similarity(eta, relabel_result$eta)
eta <- relabel_result$eta
partition$clsize <- partition$clsize[relabel_result$old_to_new] 
partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
## initial b0
beta0 <- rnorm(p, m0, s0)
## initial beta
beta_diff <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p) beta_diff[,field] <- rnorm(n, 0, s)
# beta_diff[partition$clsize==0,] <- NA
## create alpha from beta_diff and beta_0
beta_to_alpha <- function(beta_diff, beta_0){
  alpha <- beta_diff 
  for(i in 1:nrow(beta_diff)){
    # if(!is.na(beta_diff[i,1])){
    exp_beta <- exp(beta_0 + beta_diff[i,])
    alpha[i,] <- exp_beta/(exp_beta+1)
    # }
  }
  return(alpha)
}
##
alpha <- beta_to_alpha(beta_diff, beta0)
## initial frequencies of categories, initial values
theta <- fieldfrequencies[1:p]
## sanity check
all(Mvec[1:p] == sapply(theta, length))

## compute log-likelihood associated with each cluster
partition_ll <- compute_loglikelihood_all_clusters_all_fields_cpp(partition, theta, V - 1, alpha)
icluster <- 10
l <- 2

partition_ll[icluster,l]
compute_loglikelihood_one_cluster_one_field_cpp(l - 1, icluster - 1 , partition, theta[[l]], V - 1, alpha[icluster,l])
compute_loglikelihood_all_clusters_one_field_cpp(l -1 , partition, theta[[l]], V -1, alpha)[icluster]
