#'@export
single_gibbs <- function(nmcmc, V, fieldfrequencies, hyper, precomp, verbose = TRUE){
  n <- dim(V)[1]
  p <- dim(V)[2]
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## draw initial state
  state <- rinit(n, fieldfrequencies, hyper)
  
  ## record history of certain components
  N_history <- rep(NA, nmcmc)
  theta_history <- list()
  for (field in 1:p){
    theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field])
  }
  beta_0_history <- matrix(NA, nrow = nmcmc, ncol = p)
  ksize_history <- rep(NA, nmcmc)
  
  ## sanity check
  all(Mvec[1:p] == sapply(state$theta, length))
  ## compute log-likelihood associated with each cluster
  partition_ll <- couplingdeduplication:::compute_loglikelihood_all_clusters_all_fields_cpp(state$partition, state$theta, 
                                                                                            state$logtheta, V-1, state$alpha)
  ## next iterate updates in the Gibbs sampler
  for (imcmc in 1:nmcmc){
    if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize = ', state$partition$ksize, 'N = ', state$N,  "\n")
    
    update_eta_result <- couplingdeduplication:::update_eta_cpp(state$eta-1, state$partition, partition_ll, state$theta, state$logtheta, V-1, 
                                                                state$alpha, state$N)
    state$eta <- update_eta_result$eta + 1
    partition_ll <- update_eta_result$clusterloglikelihoods
    state$partition <- update_eta_result$clustering
    ksize_history[imcmc] <- state$partition$ksize
    ## 
    if (verbose) cat("NA in partition ll make sense?", all(is.na(partition_ll[,1]) == (state$partition$clsize==0)), "\n")
    
    ## relabel 
    relabel_result <- relabel(state$eta, state$partition)
    state$eta <- relabel_result$eta
    state$partition$clsize <- state$partition$clsize[relabel_result$old_to_new]
    state$partition$clmembers <- state$partition$clmembers[relabel_result$old_to_new,]
    partition_ll <- partition_ll[relabel_result$old_to_new,]
    state$alpha <- state$alpha[relabel_result$old_to_new,]
    state$beta_diff <- state$beta_diff[relabel_result$old_to_new,]
    ## update of N
    ## truncate N to N_max
    log_N_weights <- rep(-Inf, precomp$N_max)
    possiblevalues <- (state$partition$ksize):(precomp$N_max)
    log_N_weights[possiblevalues] <- precomp$lfactorials[possiblevalues+1] -
      precomp$lfactorials[possiblevalues+1-state$partition$ksize] -
      (state$n+hyper$g) * precomp$lns[possiblevalues+1]
    max_log_N_weights <- max(log_N_weights)
    N_weights <- exp(log_N_weights - max_log_N_weights)
    state$N <- sample.int(n = N_max, size = 1, prob = N_weights)
    N_history[imcmc] <- state$N
    ## update of  beta_0
    ## using the formula for the conjugate posterior in the Normal-Normal model
    common_var <- 1/(1/hyper$s0_sq + state$n/hyper$s_sq)
    for (field in 1:p){
      # calculating beta_prime as previous_beta_0 + new differences
      beta_prime <- state$beta_0[field] + state$beta_diff[,field]
      # posterior mean in the Normal-Normal model
      mu <- common_var*(hyper$m0/hyper$s0_sq + sum(beta_prime)/hyper$s_sq)
      state$beta_0[field] <- rnorm(1, mean = mu, sd = sqrt(common_var))
    }
    beta_0_history[imcmc,] <- state$beta_0
    cl_size <- state$partition$clsize
    ## update of alpha/beta
    for (icluster in 1:state$n){
      ## for each cluster...
      if (cl_size[icluster] <= 1){
        ## this cluster is empty or contains only one element so likelihood
        ## does not depend on beta parameter; thus beta follows the prior
        beta_diff_cl <- rnorm(state$p, mean = 0, sd = hyper$s)
        exp_beta_ <- exp(state$beta_0+beta_diff_cl)
        alpha_cl <- exp_beta_/(exp_beta_+1)
        # result_list <- list(beta_diff_j_l = beta_diff_j_l, alpha_j_l = alpha_j_l, accept = NA)
        state$beta_diff[icluster,] <- beta_diff_cl
        state$alpha[icluster,] <- alpha_cl
      } else {
        ## cluster has at least 2 elements, we need the likelihood into account
        # loop over fields
        for (field in 1:state$p){
          ## current beta diff
          current_beta_diff <- state$beta_diff[icluster, field]
          ## we need to compute the cluster log-likelihood
          proposed_beta_diff <- rnorm(1, current_beta_diff, sd = precomp$proposal_sd/max(1,cl_size[icluster]))
          ## calculating the corresponding value of alpha
          exp_beta <- exp(state$beta_0[field] + proposed_beta_diff)
          proposed_alpha <- exp_beta/(exp_beta+1)
          ## random number for acceptance/rejection
          logu <- log(runif(1))
          ## log likelihood taken from partition_ll + the prior
          current_posterior <-  partition_ll[icluster,field] + dnorm(current_beta_diff, 0, hyper$s, log = TRUE)
          ## we use beta_diff hence the mean is 0, not beta_0
          proposed_logprior <- dnorm(proposed_beta_diff, 0, hyper$s, log = TRUE)
          proposed_loglik <- couplingdeduplication:::compute_loglikelihood_one_cluster_one_field_cpp(field-1, icluster-1, state$partition,
                                                                                                     state$theta[[field]], state$logtheta[[field]], V-1, proposed_alpha)
          accept <- (logu < (proposed_loglik + proposed_logprior - current_posterior))
          if (accept){
            ## update log-likelihood associated cluster
            partition_ll[icluster,field] <- proposed_loglik
            state$beta_diff[icluster,field] <- proposed_beta_diff
            state$alpha[icluster,field] <- proposed_alpha
          } else {
            ## other wise do nothing
          }
        }
      }
    }
    ## update of theta
    ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
    update_theta_result <- update_theta(state$theta, state$partition, partition_ll, state$alpha, precomp$concentration)
    state$theta <- update_theta_result$theta
    state$logtheta <- lapply(state$theta, function(x) log(x))
    for (field in 1:p){
      theta_history[[field]][imcmc, ] <- state$theta[[field]]
    }
    partition_ll <- update_theta_result$partition_ll
  }
  return(list(ksize_history = ksize_history, 
              N_history = N_history,
              beta_0_history = beta_0_history,
              theta_history = theta_history))
}