#### this function runs one iterations of updating lambda
#### and tries to couple them
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020



coupleLambda <- function(j, lambda1, lambda2, p1, p2, a1, a2, N1, N2){
  clustering1 <- init_clustering(lambda1 - 1)
  cl_log_lik1 <- compute_loglikelihood_clusters(clustering1,
                                               p1,
                                               V,
                                               dimV, 
                                               a1)
  clustering2 <- init_clustering(lambda2 - 1)
  cl_log_lik2 <- compute_loglikelihood_clusters(clustering2,
                                                p2,
                                                V,
                                                dimV, 
                                                a2)
  psampq1 <- compute_psampq(j - 1, 
                            clustering1,
                            lambda1 - 1,
                            cl_log_lik1,
                            p1,
                            V,
                            dimV, 
                            a1,
                            N1)
  psampq2 <- compute_psampq(j - 1, 
                            clustering2,
                            lambda2 - 1,
                            cl_log_lik2,
                            p2,
                            V,
                            dimV, 
                            a2,
                            N2)
  js <- couple_multinomial_alt(p = psampq1, q= psampq2, n = n)
  return(js)
}




# source("getData.R")
# source("initPara.R")
# p1 <- p
# p2 <- p
# computeQcpp(j - 1 , lambda1 - 1, p, V, cumdime - 1, a1, n, H, N1)

# print(coupleLambda(j = 1, lambda1 = lambda1 , lambda2 = lambda2, p1 = p1, p2 = p2, a1 = a1, a2 = a2, N1 =N1, N2 = N2))
