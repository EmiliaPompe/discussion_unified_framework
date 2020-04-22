rm(list = ls())
set.seed(1)
source("getData.R")
# source("initPara.R")
# source("coupleMultinomial.R")
library(Rcpp)
# sourceCpp("computeNweights.cpp")
# sourceCpp("computeQcpp.cpp")
# source("coupleLambda.R")
# source("couple_multinomial_alt.R")
# source("coupleN.R")
# source("relabel.R")

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

## compute clustering
clustering <- list()
## size of each cluster
clustering$clsize <- rep(0, n)
## number of blocks
clustering$ksize <- 0
##

sourceCpp("init_clustering.cpp")
print(lambda)
clustering <- init_clustering(lambda-1)
clustering

compute_loglikelihood_clusters(clustering,
                            p,
                            V,
                            dimV, 
                            a
)


