rm(list = ls())
set.seed(1)
source("getData.R")
library(Rcpp)
sourceCpp("init_clustering.cpp")
sourceCpp("compute_loglikelihood_clusters.cpp")

## important for indexing
V <- V - 1 
## work with smaller data
V <- V[1:10,1:4]
## numbers of possibilities for each column
dimV <- dimV[1:4]
dime=as.integer(dimV)

## define dimension of V
n=as.integer(nrow(V))
Hdim=as.integer(ncol(V))
## frequencies of categories 
p=as.double(unlist(ALPHA[1:Hdim])) 
## this is a flat vector, alternatively can store p in an array, but this is what the C file is like
## note: p is theta in the paper

## a is alpha' 
a <- matrix(0.01, nrow = n, ncol = H)
g <- 1.02
N <- 50
## lambda is eta  in the paper
lambda <- sample(n, size = n, replace = T) 
print(lambda)
clustering <- init_clustering(lambda-1)
clustering$clsize

compute_loglikelihood_clusters(clustering, p, V, dimV, a)

## if we take a cluster with only one member we can verify log likelihood computation manually
# index of a cluster of size one
icluster <- which(clustering$clsize==1)[1]
# associated loglikelihood per field
compute_loglikelihood_clusters(clustering, p, V, dimV, a)[icluster,]
# index of corresponding observation
iobs <- clustering$clmembers[icluster,1]
sapply(1:Hdim, function(i) log(as.numeric(ALPHA[[i]][V[iobs+1,i]+1])))
# note how we have to be careful with the indexing since C indexes from 0 and R from 1
