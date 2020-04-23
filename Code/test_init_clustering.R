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
a <- a[,1:Hdim]
g <- 1.02
N <- 50
## lambda is eta  in the paper
lambda <- sample(n, size = n, replace = T) 
lambda <- sample(n, size = n, replace = T) 
print(lambda)
clustering <- init_clustering(lambda-1)
clustering$clsize

ll_cluster <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)
ll_cluster

## manual computation of the log likelihood in R, to check
for (icluster in 1:n){
  cat("cluster", icluster, "of size", clustering$clsize[icluster], "\n")
  if (clustering$clsize[icluster]==0){
    print(rep(-Inf,Hdim))
  } else {
    members <- clustering$clmembers[icluster,1:clustering$clsize[icluster]]
    cat("members", members, "\n")
    cat("ll computed in Rcpp:\n", ll_cluster[icluster,], "\n")
    firstmember <- members[1]
    ## compute associated log likelihood
    ll_per_field <- as.numeric(log(sapply(1:Hdim, function(index) (ALPHA[[index]])[V[firstmember+1,index]+1])))
    ## if there are more than one member, implement recursion
    if (clustering$clsize[icluster]>1){
      for (imember in 2:clustering$clsize[icluster]){
        currentmember <- members[imember]
        for (l in 1:Hdim){
          ll_per_field[l] <- ll_per_field[l] + log(a[icluster,l] * (ALPHA[[l]])[V[currentmember+1,l]+1])
          ## and for the second term loop over other members of the cluster
          secondterm <- log((1-a[icluster,l])) + log((ALPHA[[l]])[V[currentmember+1,l]+1])
          for (iothermember in 1:imember){
            othermember <- members[iothermember]
            secondterm <- secondterm + log((1-a[icluster,l])*(V[othermember+1,l]==V[currentmember+1,l]) + a[icluster,l] * (ALPHA[[l]])[V[othermember+1,l]+1])
          }
          # cat("adds second term to ", l, "-th term:", secondterm, "\n")
          max_terms <- max(secondterm, ll_per_field[l])
          ll_per_field[l] <- max_terms + log(exp(ll_per_field[l] - max_terms) + exp(secondterm - max_terms))
        }
      }
    }  
    cat("ll computer in R:\n", ll_per_field, "\n")
  }
}


