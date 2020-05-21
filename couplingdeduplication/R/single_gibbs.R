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
      state <- update_theta(state, V, algotuning)
      for (field in 1:p){
        theta_history[[field]][imcmc, ] <- state$theta[[field]]
      }
    }
    ## update N given rest
    state <- update_N(state, V, hyper, algotuning)
    N_history[imcmc] <- state$N
  }
  return(list(ksize_history = ksize_history, N_history = N_history, theta_history = theta_history))
}
