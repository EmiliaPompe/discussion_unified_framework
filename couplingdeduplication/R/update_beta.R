#'@export 
update_beta <- function(state, hyper, V, algotuning){
  n <- dim(V)[1]
  p <- dim(V)[2]
  ## using the formula for the conjugate posterior in the Normal-Normal model
  common_var <- 1/(1/hyper$s0_sq + n/hyper$s_sq)
  for (field in 1:p){
    # calculating beta_prime as previous_beta_0 + new differences
    beta_prime <- state$beta_0[field] + state$beta_diff[,field]
    # posterior mean in the Normal-Normal model
    mu <- common_var * (hyper$m0/hyper$s0_sq + sum(beta_prime)/hyper$s_sq)
    state$beta_0[field] <- rnorm(1, mean = mu, sd = sqrt(common_var))
  }
  ## update of alpha/beta
  for (icluster in 1:n){
    ## for each cluster...
    if (state$partition$clsize[icluster] <= 1){
      ## this cluster is empty or contains only one element so likelihood
      ## does not depend on beta parameter; thus beta follows the prior
      beta_diff_cl <- rnorm(p, mean = 0, sd = hyper$s)
      exp_beta_ <- exp(state$beta_0+beta_diff_cl)
      alpha_cl <- exp_beta_/(exp_beta_+1)
      state$beta_diff[icluster,] <- beta_diff_cl
      state$alpha[icluster,] <- alpha_cl
    } else {
      ## cluster has at least 2 elements, we need to take the likelihood into account
      # loop over fields
      for (field in 1:p){
        ## current beta diff
        current_beta_diff <- state$beta_diff[icluster, field]
        ## we need to compute the cluster log-likelihood
        proposed_beta_diff <- rnorm(1, current_beta_diff, sd = algotuning$proposal_sd/max(1,state$partition$cl_size[icluster]))
        ## calculating the corresponding value of alpha
        exp_beta <- exp(state$beta_0[field] + proposed_beta_diff)
        proposed_alpha <- exp_beta/(exp_beta+1)
        ## random number for acceptance/rejection
        logu <- log(runif(1))
        ## log likelihood taken from partition_ll + the prior
        current_posterior <-  state$partition_ll[icluster, field] + dnorm(current_beta_diff, 0, hyper$s, log = TRUE)
        ## we use beta_diff hence the mean is 0, not beta_0
        proposed_logprior <- dnorm(proposed_beta_diff, 0, hyper$s, log = TRUE)
        proposed_loglik <- couplingdeduplication:::compute_loglikelihood_one_cluster_one_field_cpp(field-1, icluster-1, state$partition,
                                                                                                   state$theta[[field]], state$logtheta[[field]], V-1, proposed_alpha)
        accept <- (logu < (proposed_loglik + proposed_logprior - current_posterior))
        if (accept){
          ## update log-likelihood associated cluster
          state$partition_ll[icluster,field] <- proposed_loglik
          state$beta_diff[icluster,field] <- proposed_beta_diff
          state$alpha[icluster,field] <- proposed_alpha
        } else {
          ## other wise do nothing
        }
      }
    }
  }
  return(state)
}

#'@export 
coupled_update_beta <- function(state1, state2, hyper, V, algotuning){
  n <- dim(V)[1]
  p <- dim(V)[2]
  ## using the formula for the conjugate posterior in the Normal-Normal model
  common_var <- 1/(1/hyper$s0_sq + n/hyper$s_sq)
  for (field in 1:p){
    # calculating beta_prime as previous_beta_0 + new differences
    beta_prime1 <- state1$beta_0[field] + state1$beta_diff[,field]
    beta_prime2 <- state2$beta_0[field] + state2$beta_diff[,field]
    # posterior mean in the Normal-Normal model
    mu1 <- common_var * (hyper$m0/hyper$s0_sq + sum(beta_prime1)/hyper$s_sq)
    mu2 <- common_var * (hyper$m0/hyper$s0_sq + sum(beta_prime2)/hyper$s_sq)
    beta_0s <- reflectionmaxcoupling(mu1, mu2, sqrt(common_var))
    state1$beta_0[field] <- beta_0s$xy[1]
    state2$beta_0[field] <- beta_0s$xy[2]
  }
  ## update of alpha/beta
  for (icluster in 1:n){
    ## three cases
    ## first case: in both chains, cluster size <= 1
    if ((state1$partition$clsize[icluster] <= 1) || (state2$partition$clsize[icluster] <= 1)){
      ## this cluster is empty or contains only one element so likelihood
      ## does not depend on beta parameter; thus beta follows the prior
      ## common random number strategy
      beta_diff_cl <- rnorm(p, mean = 0, sd = hyper$s)
      exp_beta_1 <- exp(state1$beta_0+beta_diff_cl)
      exp_beta_2 <- exp(state2$beta_0+beta_diff_cl)
      alpha_cl1 <- exp_beta_1/(exp_beta_1+1)
      alpha_cl2 <- exp_beta_2/(exp_beta_2+1)
      state1$beta_diff[icluster,] <- beta_diff_cl
      state1$alpha[icluster,] <- alpha_cl1
      state2$beta_diff[icluster,] <- beta_diff_cl
      state2$alpha[icluster,] <- alpha_cl2
    }
    ## second case: in both chains, cluster size >= 1
    if ((state1$partition$clsize[icluster] > 1) || (state2$partition$clsize[icluster] > 1)){
      ## clusters have at least 2 elements, we need to take the likelihood into account
      # loop over fields
      for (field in 1:p){
        ## current beta diff
        current_beta_diff1 <- state1$beta_diff[icluster, field]
        current_beta_diff2 <- state2$beta_diff[icluster, field]
        ## we need to compute the cluster log-likelihood
        proposals_ <- maxcoupling_normal(current_beta_diff1, current_beta_diff2, algotuning$proposal_sd/max(1,state1$partition$cl_size[icluster]),
                                         algotuning$proposal_sd/max(1,state2$partition$cl_size[icluster]))
        proposed_beta_diff1 <- proposals_$xy[1]
        proposed_beta_diff2 <- proposals_$xy[2]
        ## calculating the corresponding value of alpha
        exp_beta1 <- exp(state1$beta_0[field] + proposed_beta_diff1)
        exp_beta2 <- exp(state2$beta_0[field] + proposed_beta_diff2)
        proposed_alpha1 <- exp_beta1/(exp_beta1+1)
        proposed_alpha2 <- exp_beta2/(exp_beta2+1)
        ## random number for acceptance/rejection
        logu <- log(runif(1))
        ## log likelihood taken from partition_ll + the prior
        current_posterior1 <-  state1$partition_ll[icluster,field] + dnorm(current_beta_diff1, 0, hyper$s, log = TRUE)
        current_posterior2 <-  state2$partition_ll[icluster,field] + dnorm(current_beta_diff2, 0, hyper$s, log = TRUE)
        ## we use beta_diff hence the mean is 0, not beta_0
        proposed_logprior1 <- dnorm(proposed_beta_diff1, 0, hyper$s, log = TRUE)
        proposed_logprior2 <- dnorm(proposed_beta_diff2, 0, hyper$s, log = TRUE)
        proposed_loglik1 <- couplingdeduplication:::compute_loglikelihood_one_cluster_one_field_cpp(field-1, icluster-1, state1$partition,
                                                                                                   state1$theta[[field]], state1$logtheta[[field]], V-1, proposed_alpha1)
        proposed_loglik2 <- couplingdeduplication:::compute_loglikelihood_one_cluster_one_field_cpp(field-1, icluster-1, state2$partition,
                                                                                                   state2$theta[[field]], state2$logtheta[[field]], V-1, proposed_alpha2)
        accept1 <- (logu < (proposed_loglik1 + proposed_logprior1 - current_posterior1))
        accept2 <- (logu < (proposed_loglik2 + proposed_logprior2 - current_posterior2))
        if (accept1){
          ## update log-likelihood associated cluster
          state1$partition_ll[icluster,field] <- proposed_loglik1
          state1$beta_diff[icluster,field] <- proposed_beta_diff1
          state1$alpha[icluster,field] <- proposed_alpha1
        }
        if (accept2){
          ## update log-likelihood associated cluster
          state2$partition_ll[icluster,field] <- proposed_loglik2
          state2$beta_diff[icluster,field] <- proposed_beta_diff2
          state2$alpha[icluster,field] <- proposed_alpha2
        }
      }
    }
    ## third case: one chain has cluster size <= 1 and the other > 1
    if ((state1$partition$clsize[icluster] <= 1) || (state2$partition$clsize[icluster] > 1)){
      ## TODO       
    }
    if ((state1$partition$clsize[icluster] > 1) || (state2$partition$clsize[icluster] <= 1)){
      ## TODO 
    }
    
  }
  return(list(state1 = state1, state2 = state2))
}




