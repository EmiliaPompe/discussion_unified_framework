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
cumdime <- 1 +  c(0,  cumsum(dime)) 
## define dimension of V
n=as.integer(nrow(V))
Hdim=as.integer(ncol(V))
## frequencies of categories 
p=as.double(unlist(ALPHA[1:Hdim])) 
## this is a flat vector, alternatively can store p in an array, but this is what the C file is like
## note: p is theta in the paper

## a is alpha' 
a <- matrix(0.01, nrow = n, ncol = H)
a <- a[,1:Hdim]
g <- 1.02
N <- 50
## lambda is eta  in the paper
lambda <- sample(n, size = n, replace = T) 
lambda <- sample(n, size = n, replace = T) 
print(lambda)
clustering <- init_clustering(lambda-1)
clustering$clsize
ll_clusters <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)

sourceCpp("update_clustering.cpp")
j <- 1
## moves the 1st record from 5th cluster  to 4th cluster
update_clustering(clustering, j - 1, lambda[j] - 1, 3)


sourceCpp("compute_include_exclude_clustering_likelihood.cpp")
in_ex_result <- compute_include_exclude(j - 1, clustering, lambda - 1, ll_clusters, p, V, dimV, a, N)
inc <- in_ex_result$logprob_include
exc <- in_ex_result$logprob_exclude
print(inc)
print(exc)
print(all(inc <= exc))

print(in_ex_result$logprob_newcluster)
### check 
logweights <- rowSums(inc - exc) + in_ex_result$logprob_newcluster
maxlogweights <- max(logweights)
logweights <- logweights - maxlogweights
weights <- exp(logweights) / sum(exp(logweights))
weights
in_ex_result$psampq
