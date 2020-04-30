#### this function runs one iterations of updating lambda
#### and tries to couple them
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020

# coupleLambda <- function(j, lambda1, lambda2, p1, p2, a1, a2, N1, N2){
#   clustering1 <- init_clustering(lambda1 - 1)
#   cl_log_lik1 <- compute_loglikelihood_clusters(clustering1,
#                                                p1,
#                                                V,
#                                                dimV, 
#                                                a1)
#   clustering2 <- init_clustering(lambda2 - 1)
#   cl_log_lik2 <- compute_loglikelihood_clusters(clustering2,
#                                                 p2,
#                                                 V,
#                                                 dimV, 
#                                                 a2)
#   psampq1 <- compute_psampq(j - 1, 
#                             clustering1,
#                             lambda1 - 1,
#                             cl_log_lik1,
#                             p1,
#                             V,
#                             dimV, 
#                             a1,
#                             N1)
#   psampq2 <- compute_psampq(j - 1, 
#                             clustering2,
#                             lambda2 - 1,
#                             cl_log_lik2,
#                             p2,
#                             V,
#                             dimV, 
#                             a2,
#                             N2)
#   js <- couple_multinomial_alt(p = psampq1, q= psampq2, n = n)
#   return(js)
# }



couple_eta <- function(j, eta1, eta2, partition1, partition2, cl_log_lik1, cl_log_lik2, theta1, theta2, alpha1, alpha2, N1, N2){
  oldlabel1 <- eta1[j]
  oldlabel2 <- eta2[j]
  # compute everything: probability of including j excluding j from each cluster
  # full conditional of record j belong to each cluster
  res_psampq1 <- compute_include_exclude(j - 1, partition1, eta1 - 1, cl_log_lik1, theta1, V, dimV, alpha1, N1)
  res_psampq2 <- compute_include_exclude(j - 1, partition2, eta2 - 1, cl_log_lik2, theta2, V, dimV, alpha2, N2)
  psampq1 <- res_psampq1$psampq
  psampq2 <- res_psampq2$psampq
  newlabels <- couple_multinomial_alt(p = psampq1, q= psampq2, n = n)
  newlabel1 <- newlabels[1]
  newlabel2 <- newlabels[2]
  # update partitons 
  partition1 <- update_clustering(partition1, j - 1, oldlabel1 - 1, newlabel1 - 1)
  partition2 <- update_clustering(partition2, j - 1, oldlabel2 - 1, newlabel2 - 1)
  # update lambdas (must be down after update partitions, because update_clustering takes label from old lambda as input)
  eta1[j] <- newlabel1
  eta2[j] <- newlabel2
  # update cluster log likelihoods
  new_cl_log_lik1 <- cl_log_lik1
  new_cl_log_lik1[oldlabel1, ] <- res_psampq1$logprob_exclude[oldlabel1 , ]
  new_cl_log_lik1[newlabel1, ] <- res_psampq1$logprob_include[newlabel1 , ]
  new_cl_log_lik2 <- cl_log_lik2
  new_cl_log_lik2[oldlabel2, ] <- res_psampq2$logprob_exclude[oldlabel2 , ]
  new_cl_log_lik2[newlabel2, ] <- res_psampq2$logprob_include[newlabel2 , ]
  # new_cl_log_lik1 <- compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)
  # new_cl_log_lik2 <- compute_loglikelihood_clusters(partition2, p2, V, dimV, a2)
  return(list( eta1 = eta1, eta2 = eta2, 
               partition1 = partition1, partition2 = partition2, 
               cl_log_lik1 = new_cl_log_lik1, cl_log_lik2 = new_cl_log_lik2))
}

# source("getData.R")
# source("initPara.R")
# p1 <- p
# p2 <- p
# computeQcpp(j - 1 , lambda1 - 1, p, V, cumdime - 1, a1, n, H, N1)

# print(coupleLambda(j = 1, lambda1 = lambda1 , lambda2 = lambda2, p1 = p1, p2 = p2, a1 = a1, a2 = a2, N1 =N1, N2 = N2))
