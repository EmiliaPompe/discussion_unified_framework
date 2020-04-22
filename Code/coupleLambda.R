#### this function runs one iterations of updating lambda
#### and tries to couple them
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020

# library(Rcpp)
# sourceCpp("computeQcpp.cpp")
# source("coupleMultinomial.R")

coupleLambda <- function(j, lambda1, lambda2, p1, p2, a1, a2, N1, N2){
  psampq1 <- computeQcpp(j - 1, lambda1 - 1, p1, V, cumdime - 1, a1, n, H , N1)
  psampq2 <- computeQcpp(j - 1, lambda2 - 1, p2, V, cumdime - 1, a2, n, H , N2)
  js <- couple_multinomial_alt(p = psampq1, q= psampq2, n = n)
  return(js)
}




