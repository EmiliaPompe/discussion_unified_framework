#'@export
single_gibbs <- function(nmcmc, V, fieldfrequencies, hyper, precomp, update.theta = FALSE, verbose = TRUE){
  ## dimension of the data
  ## number of rows
  n <- dim(V)[1]
  ## number of fields
  p <- dim(V)[2]
  ## number of possibilities in each field
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## draw initial state
  state <- rinit(n, fieldfrequencies, hyper)
  ## over-write alpha, beta_0 and beta_diff 
  state$alpha <- matrix(0.01, nrow = n, ncol = p)
  state$beta_0 <- NA
  state$beta_diff <- matrix(log(0.01/0.99), nrow = n, ncol = p)
  ## record history of some of the components
  N_history <- rep(NA, nmcmc)
  ksize_history <- rep(NA, nmcmc)
  theta_history <- list()
  for (field in 1:p){ theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field]) }
  ## compute log-likelihood associated with each cluster
  partition_ll <- couplingdeduplication:::compute_loglikelihood_all_clusters_all_fields_cpp(state$partition, state$theta, 
                                                                                            state$logtheta, V-1, state$alpha)
  ## next iterate updates in the Gibbs sampler
  for (imcmc in 1:nmcmc){
    if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize = ', 
                                          state$partition$ksize, 'N = ', state$N,  "\n")
    ## update eta given the rest
    update_eta_result <- couplingdeduplication:::update_eta(state$eta-1, state$partition, partition_ll, 
                                                            state$theta, state$logtheta, V-1, 
                                                            state$alpha, state$N, 0.1)
    state$eta <- update_eta_result$eta + 1
    partition_ll <- update_eta_result$clusterloglikelihoods
    state$partition <- update_eta_result$clustering
    ksize_history[imcmc] <- state$partition$ksize
    ## relabel eta 
    relabel_result <- relabel(state$eta, state$partition)
    state$eta <- relabel_result$eta
    state$partition$clsize <- state$partition$clsize[relabel_result$old_to_new]
    state$partition$clmembers <- state$partition$clmembers[relabel_result$old_to_new,]
    partition_ll <- partition_ll[relabel_result$old_to_new,]
    state$alpha <- state$alpha[relabel_result$old_to_new,]
    state$beta_diff <- state$beta_diff[relabel_result$old_to_new,]
    ## update of theta
    if (update.theta){
      ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
      ## for each field
      for (field in 1:p){
        concentration <- precomp$concentration
        ## current state
        x <- state$theta[[field]] 
        cl_log_lik <- partition_ll[,field] 
        ## dirichlet proposal
        x_propose <- rgamma(n = length(x), shape = concentration*x+1, rate = 1)
        x_propose <- x_propose/sum(x_propose)
        logx_propose <- log(x_propose)
        cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
        ## transition ratio 
        lratio <- sum((concentration*x_propose) * log(x) - (concentration*x) * log(x_propose)) + sum(lgamma(concentration*x+1) - lgamma(concentration*x_propose+1))
        ## log accept probability
        laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
        if (log(runif(1)) < laccept){
          state$theta[[field]] <- x_propose
          state$logtheta[[field]] <- log(x_propose)
          partition_ll[,field] <- cl_log_lik_new
        }
        theta_history[[field]][imcmc, ] <- state$theta[[field]]
      }
    }
    ## update N given rest
    ## using truncation to N_max
    log_N_weights <- rep(-Inf, precomp$N_max)
    possiblevalues <- (state$partition$ksize):(precomp$N_max)
    log_N_weights[possiblevalues] <- precomp$lfactorials[possiblevalues+1] -
      precomp$lfactorials[possiblevalues+1-state$partition$ksize] -
      (state$n+hyper$g) * precomp$lns[possiblevalues+1]
    max_log_N_weights <- max(log_N_weights)
    N_weights <- exp(log_N_weights - max_log_N_weights)
    state$N <- sample.int(n = N_max, size = 1, prob = N_weights)
    N_history[imcmc] <- state$N
  }
  return(list(ksize_history = ksize_history, N_history = N_history, theta_history = theta_history))
}
