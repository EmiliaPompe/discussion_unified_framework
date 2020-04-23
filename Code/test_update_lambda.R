rm(list = ls())
set.seed(1)
source("getData.R")
sourceCpp("init_clustering.cpp")
sourceCpp("compute_loglikelihood_clusters.cpp")
sourceCpp("update_lambda.cpp")
library(Rcpp)
## important for indexing
V <- V - 1 
## work with smaller data
V <- V[1:10,1:4]
## numbers of possibilities for each column
dimV <- dimV[1:4]
dime=as.integer(dimV)
##
## define dimension of V
n=as.integer(nrow(V))
Hdim=as.integer(ncol(V))
## distortion parameters
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

##

print(lambda)
clustering <- init_clustering(lambda-1)
ll_clusters <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)
## now update of lambda
print(lambda)
print(clustering$clsize)
print(ll_clusters)
##
res_update <- update_lambda(lambda-1, clustering, ll_clusters, p, V, dimV, a, N)
lambda <- res_update$lambda+1
clustering <- res_update$clustering
ll_clusters <- res_update$clusterloglikelihoods

niter <- 10000
for (iter in 1:niter){
  res_update <- update_lambda(lambda-1, clustering, ll_clusters, p, V, dimV, a, N)
  lambda <- res_update$lambda+1
  clustering <- res_update$clustering
  ll_clusters <- res_update$clusterloglikelihoods
}


# # ## manual computation for cluster with 9th obs
# # log(sapply(1:Hdim, function(index) (ALPHA[[index]])[V[9,index]+1]))
# ## manual computation for cluster with 9th obs and then 4th observation
# ll_per_field <- as.numeric(log(sapply(1:Hdim, function(index) (ALPHA[[index]])[V[9,index]+1])))
# ## adds 4th obs
# icluster <- 1
# for (l in 1:Hdim){
#   ll_per_field[l] <- log(a[icluster,l] * (ALPHA[[l]])[V[4,l]+1]) + ll_per_field[l]
#   secondterm <- log((1-a[icluster,l]) * (ALPHA[[l]])[V[4,l]+1]) + log((1-a[icluster,l])*(V[4,l]==V[9,l]) + a[icluster,l] * (ALPHA[[l]])[V[9,l]+1])
#   ll_per_field[l] <- log(exp(ll_per_field[l]) + exp(secondterm))
# }
# ll_per_field
# 
# ## manual computation for cluster with 4th and 9th obs
# ll_per_field <- as.numeric(log(sapply(1:Hdim, function(index) (ALPHA[[index]])[V[4,index]+1])))
# ## adds 9th obs
# icluster <- 1
# for (l in 1:Hdim){
#   ll_per_field[l] <- log(a[icluster,l] * (ALPHA[[l]])[V[9,l]+1]) + ll_per_field[l]
#   secondterm <- log((1-a[icluster,l])) + log((ALPHA[[l]])[V[9,l]+1]) + log((1-a[icluster,l])*(V[9,l]==V[4,l]) + a[icluster,l] * (ALPHA[[l]])[V[4,l]+1])
#   ll_per_field[l] <- log(exp(ll_per_field[l]) + exp(secondterm))
# }
# ll_per_field
# ## compared to ...
# ll_clusters[1,]
# compute_loglikelihood_clusters(clustering, p, V, dimV, a)[1,]
# 
# cluster_alt <- clustering
# cluster_alt$clmembers[1,1:2] <- rev(cluster_alt$clmembers[1,1:2])
# compute_loglikelihood_clusters(cluster_alt, p, V, dimV, a)[1,]
# 
# 
# 
# ##
# cumdime <- rep(0, Hdim) 
# for (l in 2:Hdim){
#   cumdime[l] = cumdime[l-1] + dimV[l-1];
# }
# 
# ## field (R notation)
# l <- sample.int(Hdim, 1)
# ## value of observation (C notation)
# V_ <- sample.int(dimV[l], 1)-1
# ALPHA[[l]][V_+1]
# p[cumdime[l] + V_+1]
# 
# summary(as.numeric(res_update$clusterloglikelihoods) - as.numeric(compute_loglikelihood_clusters(res_update$clustering, p, V, dimV, a)))
# 
# 
# clustering$clsize
# clustering$clmembers[4,]
# compute_loglikelihood_clusters(res_update$clustering, p, V, dimV, a)[4,]
# as.numeric(log(sapply(1:Hdim, function(index) (ALPHA[[index]])[V[2,index]+1])))
