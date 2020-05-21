#'@title update_theta
#'@description update 'theta' component of the chain
#' using a mixture of random walk and independent MH 
#'@return a state with updated 'theta' and updated loglikelihood
#'@export
update_theta <- function(state, V, algotuning){
  p <- dim(V)[2]
  ## 
  u_ <- runif(1)
  if (u_ < algotuning$theta_update_rw){
    ## random walk proposal
    for (field in 1:p){
      conc <- algotuning$theta_update_rwconcentration
      ## current state
      x <- state$theta[[field]] 
      cl_log_lik <- state$partition_ll[,field] 
      ## dirichlet proposal
      x_propose <- rgamma(n = length(x), shape = conc*x+1, rate = 1)
      x_propose <- x_propose/sum(x_propose)
      logx_propose <- log(x_propose)
      cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
      ## transition ratio 
      lratio <- sum((conc*x_propose) * log(x) - (conc*x) * log(x_propose)) + sum(lgamma(conc*x+1) - lgamma(conc*x_propose+1))
      ## log accept probability
      laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
      if (log(runif(1)) < laccept){
        state$theta[[field]] <- x_propose
        state$logtheta[[field]] <- log(x_propose)
        state$partition_ll[,field] <- cl_log_lik_new
      }
    }
  } else {
    ## independent proposal
    for (field in 1:p){
      alpha_proposal <- algotuning$fieldfrequencies[[field]] * state$partition$ksize * 0.9 + 1
      ## sample from Dirichlet distribution 
      x_propose <- gtools::rdirichlet(1, alpha = alpha_proposal)
      logx_propose <- log(x_propose)
      cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
      ## current state
      x <- state$theta[[field]] 
      cl_log_lik <- state$partition_ll[,field] 
      ## transition ratio 
      lratio <- log(gtools::ddirichlet(x = x, alpha = alpha_proposal)) - log(gtools::ddirichlet(x = x_propose, alpha = alpha_proposal)) 
      ## log accept probability
      laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
      if (log(runif(1)) < laccept){
        state$theta[[field]] <- x_propose
        state$logtheta[[field]] <- log(x_propose)
        state$partition_ll[,field] <- cl_log_lik_new
      }
    }
  }
  return(state)
}


#'@title coupled_update_theta
#'@description coupled update of 'theta' component of two chains
#' using a mixture of random walk and independent MH, with common random numbers and maximal coupling ideas
#'@return a list with 'state1' and 'state2' 
#'@export
coupled_update_theta <- function(state1, state2, V, algotuning){
  p <- dim(V)[2]
  ## common random number to decide what type of update to perform
  u_ <- runif(1)
  all_theta_equal <- TRUE
  if (u_ < algotuning$theta_update_rw){
    ## random walk proposal
    for (field in 1:p){
      conc <- algotuning$theta_update_rwconcentration
      ## current state
      x1 <- state1$theta[[field]] 
      cl_log_lik1 <- state1$partition_ll[,field] 
      x2 <- state2$theta[[field]] 
      cl_log_lik2 <- state2$partition_ll[,field] 
      ## dirichlet proposal
      x_proposals <- maxcoupling_dirichlet(conc*x1+1, conc*x2+1)
      x_propose1 <- x_proposals$xy[,1]
      x_propose2 <- x_proposals$xy[,2]
      logx_propose1 <- log(x_propose1)
      logx_propose2 <- log(x_propose2)
      cl_log_lik_new1 <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state1$partition, x_propose1, logx_propose1, V - 1, state1$alpha)
      cl_log_lik_new2 <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state2$partition, x_propose2, logx_propose2, V - 1, state2$alpha)
      ## transition ratio 
      lratio1 <- sum((conc*x_propose1) * log(x1) - (conc*x1) * log(x_propose1)) + sum(lgamma(conc*x1+1) - lgamma(conc*x_propose1+1))
      lratio2 <- sum((conc*x_propose2) * log(x2) - (conc*x2) * log(x_propose2)) + sum(lgamma(conc*x2+1) - lgamma(conc*x_propose2+1))
      ## log accept probability
      laccept1 <- lratio1 + sum(cl_log_lik_new1[which(state1$partition$clsize != 0)]) - sum(cl_log_lik1[which(state1$partition$clsize != 0)])
      laccept2 <- lratio2 + sum(cl_log_lik_new2[which(state2$partition$clsize != 0)]) - sum(cl_log_lik2[which(state2$partition$clsize != 0)])
      uacceptrw <- runif(1)
      if (log(uacceptrw) < laccept1){
        state1$theta[[field]] <- x_propose1
        state1$logtheta[[field]] <- log(x_propose1)
        state1$partition_ll[,field] <- cl_log_lik_new1
      }
      if (log(uacceptrw) < laccept2){
        state2$theta[[field]] <- x_propose2
        state2$logtheta[[field]] <- log(x_propose2)
        state2$partition_ll[,field] <- cl_log_lik_new2
      }
      all_theta_equal <- all_theta_equal && (x_proposals$identical && (log(uacceptrw) < laccept1) && (log(uacceptrw) < laccept2))
    }
  } else {
    ## independent proposal
    for (field in 1:p){
      alpha_proposal <- algotuning$fieldfrequencies[[field]] * min(state1$partition$ksize, state2$partition$ksize) * 0.9 + 1
      ## sample from Dirichlet distribution 
      x_propose <- gtools::rdirichlet(1, alpha = alpha_proposal)
      logx_propose <- log(x_propose)
      cl_log_lik_new1 <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state1$partition, x_propose, logx_propose, V - 1, state1$alpha)
      cl_log_lik_new2 <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state2$partition, x_propose, logx_propose, V - 1, state2$alpha)
      ## current state
      x1 <- state1$theta[[field]] 
      cl_log_lik1 <- state1$partition_ll[,field] 
      x2 <- state2$theta[[field]] 
      cl_log_lik2 <- state2$partition_ll[,field] 
      ## transition ratio 
      lratio1 <- log(gtools::ddirichlet(x = x1, alpha = alpha_proposal)) - log(gtools::ddirichlet(x = x_propose, alpha = alpha_proposal)) 
      lratio2 <- log(gtools::ddirichlet(x = x2, alpha = alpha_proposal)) - log(gtools::ddirichlet(x = x_propose, alpha = alpha_proposal)) 
      ## log accept probability
      laccept1 <- lratio1 + sum(cl_log_lik_new1[which(state1$partition$clsize != 0)]) - sum(cl_log_lik1[which(state1$partition$clsize != 0)])
      laccept2 <- lratio2 + sum(cl_log_lik_new2[which(state2$partition$clsize != 0)]) - sum(cl_log_lik2[which(state2$partition$clsize != 0)])
      uacceptindep <- runif(1)
      if (log(ruacceptindep) < laccept1){
        state1$theta[[field]] <- x_propose
        state1$logtheta[[field]] <- log(x_propose)
        state1$partition_ll[,field] <- cl_log_lik_new1
      }
      if (log(ruacceptindep) < laccept2){
        state2$theta[[field]] <- x_propose
        state2$logtheta[[field]] <- log(x_propose)
        state2$partition_ll[,field] <- cl_log_lik_new2
      }
      all_theta_equal <- all_theta_equal && (log(ruacceptindep) < laccept1) && (log(ruacceptindep) < laccept2)
    }
  }
  return(list(state1 = state1, state2 = state2, identical = all_theta_equal))
}

