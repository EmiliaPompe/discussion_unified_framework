rm(list = ls())
source("getData.R")
source("initPara.R")
library(Rcpp)
sourceCpp("computeQcpp.cpp")

j <- 1
psampq <- computeQcpp(j - 1, lambda - 1, p, V, cumdime - 1, a, n, H , N1)

a1 <- a 
lambda1 <- lambda 
N1 <- N1 
p1 <- p

a2 <- matrix(0.02, nrow = n, ncol = H)
lambda2 <- sample(n, size = n, replace = TRUE)
p2 <- p
N2 <- N1 + 10 


print( coupleLambda(j = j, lambda1 = lambda1 , lambda2 = lambda2, p1 = p1, p2 = p2, a1 = a1, a2 = a2, N1 =N1, N2 = N2))




