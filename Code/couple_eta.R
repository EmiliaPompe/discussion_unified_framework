#### this function runs one iterations of updating eta_j for two chains
#### such that each update marginally follows the full conditional distribution 
#### and that the two chains are coupled
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 11-May 2020

couple_eta <- function(j, eta1, eta2, partition1, partition2, partition_ll1, partition_ll2, theta1, theta2, alpha1, alpha2, N1, N2){
  oldlabel1 <- eta1[j]
  oldlabel2 <- eta2[j]
  # compute everything: probability of including j excluding j from each cluster
  # full conditional of record j belong to each cluster
  res_psampq1 <- compute_include_exclude_cpp(j - 1, partition1, eta1 - 1, partition_ll1, theta1, V - 1, alpha1, N1)
  res_psampq2 <- compute_include_exclude_cpp(j - 1, partition2, eta2 - 1, partition_ll2, theta2, V - 1, alpha2, N2)
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
  new_cl_log_lik1 <- partition_ll1
  new_cl_log_lik1[oldlabel1, ] <- res_psampq1$logprob_exclude[oldlabel1 , ]
  new_cl_log_lik1[newlabel1, ] <- res_psampq1$logprob_include[newlabel1 , ]
  new_cl_log_lik2 <- partition_ll2
  new_cl_log_lik2[oldlabel2, ] <- res_psampq2$logprob_exclude[oldlabel2 , ]
  new_cl_log_lik2[newlabel2, ] <- res_psampq2$logprob_include[newlabel2 , ]
  return(list( eta1 = eta1, eta2 = eta2, 
               partition1 = partition1, partition2 = partition2, 
               partition_ll1 = new_cl_log_lik1, partition_ll2 = new_cl_log_lik2))
}

