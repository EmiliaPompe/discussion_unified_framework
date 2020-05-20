## update eta given the rest, and perform "relabeling"
#'@export
update_eta_relabel <- function(state, V, algotuning){
  ## update eta using C++ function
  update_eta_result <- couplingdeduplication:::update_eta(eta = state$eta-1, clustering = state$partition, previous_clusterloglikelihoods = state$partition_ll, 
                                                          theta = state$theta, logtheta = state$logtheta, V = V-1, 
                                                          alpha = state$alpha, N = state$N, updateprobability = algotuning$eta_update_prob)
  state$eta <- update_eta_result$eta + 1
  state$partition_ll <- update_eta_result$clusterloglikelihoods
  state$partition <- update_eta_result$clustering
  ## relabel eta 
  relabel_result <- relabel(state$eta, state$partition)
  state$eta <- relabel_result$eta
  state$partition$clsize <- state$partition$clsize[relabel_result$old_to_new]
  state$partition$clmembers <- state$partition$clmembers[relabel_result$old_to_new,]
  state$partition_ll <- state$partition_ll[relabel_result$old_to_new,]
  state$alpha <- state$alpha[relabel_result$old_to_new,]
  state$beta_diff <- state$beta_diff[relabel_result$old_to_new,]
  return(state)
}


#'@export
coupled_update_eta_relabel <- function(state1, state2, V, algotuning){
  ## draw new etas
  update_eta_result <- couplingdeduplication:::coupled_update_eta(state1$eta-1, state2$eta-1, state1$partition, state2$partition, 
                                                                  state1$partition_ll, state2$partition_ll, state1$theta, state2$theta, 
                                                                  state1$logtheta, state2$logtheta, V-1, 
                                                                  state1$alpha, state2$alpha, state1$N, state2$N, algotuning$eta_update_prob)
  state1$eta <- update_eta_result$eta1 + 1
  state1$partition_ll <- update_eta_result$clusterloglikelihoods1
  state1$partition <- update_eta_result$clustering1
  state2$eta <- update_eta_result$eta2 + 1
  state2$partition_ll <- update_eta_result$clusterloglikelihoods2
  state2$partition <- update_eta_result$clustering2
  ## relabel 
  relabel_result1 <- relabel(state1$eta, state1$partition)
  state1$eta <- relabel_result1$eta
  state1$partition$clsize <- state1$partition$clsize[relabel_result1$old_to_new]
  state1$partition$clmembers <- state1$partition$clmembers[relabel_result1$old_to_new,]
  state1$partition_ll <- state1$partition_ll[relabel_result1$old_to_new,]
  state1$alpha <- state1$alpha[relabel_result1$old_to_new,]
  state1$beta_diff <- state1$beta_diff[relabel_result1$old_to_new,]
  #
  relabel_result2 <- relabel(state2$eta, state2$partition)
  state2$eta <- relabel_result2$eta
  state2$partition$clsize <- state2$partition$clsize[relabel_result2$old_to_new]
  state2$partition$clmembers <- state2$partition$clmembers[relabel_result2$old_to_new,]
  state2$partition_ll <- state2$partition_ll[relabel_result2$old_to_new,]
  state2$alpha <- state2$alpha[relabel_result2$old_to_new,]
  state2$beta_diff <- state2$beta_diff[relabel_result2$old_to_new,]
  return(list(state1 = state1, state2 = state2))
}
