## Script to try the maximal coupling of Dirichlet distributions
## this is work in progress
rm(list = ls())
set.seed(2)
library(couplingdeduplication)

size_ <- 3
alpha1 <- rexp(size_) 
alpha2 <- rexp(size_) 
alphatilde1 <- alpha1 / sum(alpha1)
alphatilde2 <- alpha2 / sum(alpha2)

nrep <- 1e4
samples1 <- matrix(NA, nrow = nrep, ncol = size_)
samples2 <- matrix(NA, nrow = nrep, ncol = size_)
for (irep in 1:nrep){
  results <- maxcoupling_dirichlet(alpha1, alpha2)
  samples1[irep,] <- results$xy[,1]
  samples2[irep,] <- results$xy[,2]
}

colMeans(samples1)
alphatilde1
colMeans(samples2)
alphatilde2

apply(samples1, 2, function(v) var(v))
alphatilde1 * (1 - alphatilde1) / (1+sum(alpha1))

apply(samples2, 2, function(v) var(v))
alphatilde2 * (1 - alphatilde2) / (1+sum(alpha2))

