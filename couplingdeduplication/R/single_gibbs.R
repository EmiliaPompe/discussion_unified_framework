#'@title single_gibbs
#'@description runs a Gibbs sampler with updates of eta, N, and optionally of theta
#'@return a list with  'ksize_history', 'N_history' and 'theta_history'
#'@export
single_gibbs <- function(nmcmc, V, fieldfrequencies, hyper, algotuning, update.theta = FALSE){
  ## dimension of the data
  ## number of rows
  n <- dim(V)[1]
  ## number of fields
  p <- dim(V)[2]
  ## number of possibilities in each field
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## draw initial state
  state <- rinit(n, fieldfrequencies, hyper, V)
  ## record history of some of the components
  N_history <- rep(NA, nmcmc)
  ksize_history <- rep(NA, nmcmc)
  theta_history <- list()
  for (field in 1:p){ theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field]) }
  ## next iterate updates in the Gibbs sampler
  for (imcmc in 1:nmcmc){
    if (algotuning$verbose && (imcmc %% 10 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize = ', 
                                          state$partition$ksize, 'N = ', state$N,  "\n")
    ## update eta given the rest
    state <- update_eta_relabel(state, V, algotuning)
    ksize_history[imcmc] <- state$partition$ksize
    ## update of theta
    if (update.theta){
      ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
      ## for each field
      for (field in 1:p){
        theta_update_conc <- algotuning$theta_update_conc
        ## current state
        x <- state$theta[[field]] 
        cl_log_lik <- state$partition_ll[,field] 
        ## dirichlet proposal
        x_propose <- rgamma(n = length(x), shape = theta_update_conc*x+1, rate = 1)
        x_propose <- x_propose/sum(x_propose)
        logx_propose <- log(x_propose)
        cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
        ## transition ratio 
        lratio <- sum((theta_update_conc*x_propose) * log(x) - (theta_update_conc*x) * log(x_propose)) + sum(lgamma(theta_update_conc*x+1) - lgamma(theta_update_conc*x_propose+1))
        ## log accept probability
        laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
        if (log(runif(1)) < laccept){
          state$theta[[field]] <- x_propose
          state$logtheta[[field]] <- log(x_propose)
          state$partition_ll[,field] <- cl_log_lik_new
        }
        theta_history[[field]][imcmc, ] <- state$theta[[field]]
      }
    }
    ## update N given rest
    state <- update_N(state, V, hyper, algotuning)
    N_history[imcmc] <- state$N
  }
  return(list(ksize_history = ksize_history, N_history = N_history, theta_history = theta_history))
}
