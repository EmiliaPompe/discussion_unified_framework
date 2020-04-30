rm(list = ls())
set.seed(1)
test <- T
source("getData.R")
source("initPara.R")
# source("coupleMultinomial.R")
library(Rcpp)
# sourceCpp("computeNweights.cpp")
# sourceCpp("computeQcpp.cpp")
# source("coupleLambda.R")
# source("couple_multinomial_alt.R")
# source("coupleN.R")
# source("relabel.R")
sourceCpp("init_clustering.cpp")
sourceCpp("compute_loglikelihood_clusters.cpp")
sourceCpp('compute_loglikelihood_clusters_field.cpp')
sourceCpp('compute_psampq.cpp')
sourceCpp("update_clustering.cpp")
sourceCpp("compute_include_exclude_clustering_likelihood.cpp")
source("couple_multinomial_alt.R")
source("coupleEta.R")
source("relabel.R")
eta1 <- relabel(eta1)$lambda
eta2 <- relabel(eta2)$lambda

eta <- sample(n, size = n, replace = T) 
print(eta)
clustering <- init_clustering(eta-1)
clustering$clsize

cl_log_lik <- compute_loglikelihood_clusters(clustering, theta, V, dimV, alpha)
j <- 3
qs <- compute_psampq(j - 1, 
                     clustering,
                     eta - 1,
                     cl_log_lik,
                     theta,
                     V,
                     dimV, 
                     alpha,
                     N1)


j <- 2

partition1 <- init_clustering(eta1 - 1)
cl_log_lik1 <- compute_loglikelihood_clusters(partition1,
                                              theta1,
                                              V,
                                              dimV, 
                                              alpha1)
partition2 <- init_clustering(eta2 - 1)
cl_log_lik2 <- compute_loglikelihood_clusters(partition2,
                                              theta2,
                                              V,
                                              dimV, 
                                              alpha2)

result <- couple_eta(j, eta1, eta2, partition1, partition2, cl_log_lik1, cl_log_lik2, theta1, theta2, alpha1, alpha2, N1, N2)

## new cluster probabitliy
result$cl_log_lik1
## alternatively
compute_loglikelihood_clusters(result$partition1, theta1, V, dimV, alpha1)
all(result$cl_log_lik1 == compute_loglikelihood_clusters(result$partition1, theta1, V, dimV, alpha1))

for(j in 1:n){
  result <- couple_eta(j, eta1, eta2, partition1, partition2, cl_log_lik1, cl_log_lik2, theta1, theta2, alpha1, alpha2, N1, N2)
  print( (result$cl_log_lik1 - compute_loglikelihood_clusters(result$partition1, theta1, V, dimV, alpha1)) )
}

