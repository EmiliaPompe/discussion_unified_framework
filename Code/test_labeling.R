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
V <- V[1:40,1:4]
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
clustering <- init_clustering(lambda-1)
ll_clusters <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)

niter <- 1000
for (iter in 1:niter){
  res_update <- update_lambda(lambda-1, clustering, ll_clusters, p, V, dimV, a, N)
  lambda <- res_update$lambda+1
  clustering <- res_update$clustering
  ll_clusters <- res_update$clusterloglikelihoods
}


## returns new labels "newlambda"
## and correspondence between old labels and new labels in "iis"
## i.e. iis[n] gives 
relabel1 <- function(lambda){
  n <- length(lambda)
  iis <- 1:n
  newlambda <- rep(NA, n)
  visited <- rep(FALSE, n)
  changed <- rep(FALSE, n)
  currentCnt <- 1
  for (j in 1:n){
    if (!visited[j]){
      iis[lambda[j]] <- currentCnt
      changed[lambda[j]] <- TRUE
      clMembers <- which(lambda == lambda[j])
      visited[clMembers] <- TRUE
      newlambda[clMembers] <- currentCnt
      currentCnt <- currentCnt + 1
    }
  }
  iis[!changed] <- c(currentCnt : n)
  return(list(lambda = newlambda, iis = iis))
}

lambda
relabelled_lambda <- relabel1(lambda)
lambda[4]
relabelled_lambda$iis[4]
relabelled_lambda$lambda[relabelled_lambda$iis[4]]

cl1 <- init_clustering(lambda-1)
cl2 <- init_clustering(relabelled_lambda$lambda-1)

## one way of checking this is the same clustering is to compare
## comembership matrix
clusteval::cluster_similarity(lambda, relabelled_lambda$lambda)


## function to perform re-labeling
## this one goes through each cluster and gives it the label
## of smallest index of member of cluster
relabel2 <- function(lambda, clustering){
  newlambda <- rep(0, length(lambda))
  ## permutation is a vector with i-th element equal to the new label given to what 
  ## was previously the i-th label 
  permutation <- rep(0, length(lambda))
  ## for each cluster...
  for (icluster in 1:length(clustering$clsize)){
    ## if cluster is not empty
    if (clustering$clsize[icluster]>0){
      ## get indices of cluster members
      cluster_members <- clustering$clmembers[icluster,1:clustering$clsize[icluster]]
      ## give as label smallest of these indices
      cluster_label <- min(cluster_members)
      newlambda[cluster_members+1] <- cluster_label + 1
      permutation[icluster] <- cluster_label + 1
    }
  }
  ## assign arbitrary labels to empty clusters
  permutation[permutation==0] <- setdiff(1:length(lambda), sort(permutation[permutation!=0]))
  # newcl <- init_clustering(newlambda-1)
  # identical(clustering$clsize[order(permutation)], newcl$clsize)
  # identical(clustering$clmembers[order(permutation),], newcl$clmembers)
  # identical(clustering$ksize, newcl$ksize)
  return(list(lambda = newlambda, clustering = list(clsize = clustering$clsize[order(permutation)], ksize = clustering$ksize,
                                                    clmembers = clustering$clmembers[order(permutation),])))
}



lambda <- sample(n, size = n, replace = T) 
clustering <- init_clustering(lambda-1)
result_ <- relabel2(lambda, clustering)

clusteval::cluster_similarity(lambda, result_$lambda)

