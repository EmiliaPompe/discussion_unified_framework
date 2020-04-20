#### this function runs one iterations of coupled update of N
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020


library(Rcpp)
sourceCpp("computeNweights.cpp")
source("coupleMultinomial.R")

N1 <- 15
N2 <- 20
stepsize <- 3
ksize1 <- 10
ksize2 <- 13
g <- 1.02
n <- 4
coupleN <- function(N1, N2, k1,k2, g, n){
  w <- computeNweights(g = g, N1 =N1, N2  = N2, stepsize = stepsize, ksize1 = ksize1,ksize2 = ksize2, n = n)
  index <- coupleMultinomial(p = w[2,],q = w[3,], n = dim(w)[2])
  return(as.integer(w[1,index]))
}

coupleN(N1, N2, k1, k2, g, n)
