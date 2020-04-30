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
source("couple_dirichlet.R")
source("coupleTheta.R")


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

count <- 0 
distance <- rep(0, 100)
for (iter in 1:100){
  couple_result <- coupleTheta(theta1, theta2, partition1, partition2, alpha1, alpha2)
  count <- count + all(couple_result$theta1 == couple_result$theta2)
  distance[iter] <- sum((couple_result$theta1 - couple_result$theta2)**2)
}
count / 100
plot(distance)

concentration <- 10**4
distance <- rep(0, 1000)
for( iter in 1:1000){
  couple_result <- coupleTheta(theta1, theta2, partition1, partition2, alpha1, alpha2)
  theta1 <- couple_result$theta1
  theta2 <- couple_result$theta2
  distance[iter] <- sum((theta1 - theta2)**2)
}
plot(distance, type = 'l')

