#### this function runs one iterations of coupled update of N
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020

# 
# library(Rcpp)
# sourceCpp("computeNweights.cpp")
# source("coupleMultinomial.R")

couple_N <- function(N1, N2, k1,k2, g, n){
  # w <- computeNweights(g = g, N1 =N1, N2  = N2, stepsize = stepsize, ksize1 = k1,ksize2 = k2, n = n)
  # index <- couple_multinomial_alt(p = w[2,],q = w[3,], n = dim(w)[2])
  # return(as.integer(w[1,index]))
  lweights1 <- sapply(1:N_max, function(xx) lfactorial(xx) - lfactorial(xx - k1) - (n+g) * log(xx) )
  lweights2 <- sapply(1:N_max, function(xx) lfactorial(xx) - lfactorial(xx - k2) - (n+g) * log(xx) )
  maxlogweights1 <- max(lweights1)
  lweights1 <- lweights1 - maxlogweights1
  w1 <- exp(lweights1) / sum(exp(lweights1))
  maxlogweights2 <- max(lweights2)
  lweights2 <- lweights2 - maxlogweights2
  w2 <- exp(lweights2) / sum(exp(lweights2))
  cp_n_res <- couple_multinomial_alt(p = w1, q = w2, n = N_max)
  return(cp_n_res)
}
