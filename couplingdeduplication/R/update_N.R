#'@title update_N
#'@description takes a state of the Markov chain and draws a new N.
#' The procedure truncates the distribution to 'N_max', computes the probabilities, and uses the 'sample.int' function of R
#'@param state a list representing the state of the Markov chain as e.g. generate by \code{\link{rinit}}
#'@param V a matrix representing the data
#'@param hyper a list containing hyperparameter 'g'
#'@param algotuning a list containing pre-computed quantities 'N_max', 'lfactorials', 'lns'
#'@return state of the Markov chain with updated 'N'
#'@export
update_N <- function(state, V, hyper, algotuning){
  ## using truncation to N_max
  log_N_weights <- rep(-Inf, algotuning$N_max)
  possiblevalues <- (state$partition$ksize):(algotuning$N_max)
  log_N_weights[possiblevalues] <- algotuning$lfactorials[possiblevalues+1] -
    algotuning$lfactorials[possiblevalues+1-state$partition$ksize] -
    (state$n+hyper$g) * algotuning$lns[possiblevalues+1]
  max_log_N_weights <- max(log_N_weights)
  N_weights <- exp(log_N_weights - max_log_N_weights)
  state$N <- sample.int(n = N_max, size = 1, prob = N_weights) 
  return(state)
}

#'@title coupled_update_N
#'@description samples from a maximal coupling of the conditional distributions of N given the other variables,
#' for the two chains given as 'state1' and 'state2'; see \code{\link{update_N}} for the other parameters.
#'@return a list with 'state1' and 'state2', the two new states of the Markov chains
#'@export
coupled_update_N <- function(state1, state2, V, hyper, algotuning){
  ## using truncation to N_max
  log_N_weights1 <- rep(-Inf, algotuning$N_max)
  log_N_weights2 <- rep(-Inf, algotuning$N_max)
  possiblevalues1 <- (state1$partition$ksize):(algotuning$N_max)
  possiblevalues2 <- (state2$partition$ksize):(algotuning$N_max)
  log_N_weights1[possiblevalues1] <- algotuning$lfactorials[possiblevalues1+1] -
    algotuning$lfactorials[possiblevalues1+1-state1$partition$ksize] -
    (state1$n+hyper$g) * algotuning$lns[possiblevalues1+1]
  log_N_weights2[possiblevalues2] <- algotuning$lfactorials[possiblevalues2+1] -
    algotuning$lfactorials[possiblevalues2+1-state2$partition$ksize] -
    (state2$n+hyper$g) * algotuning$lns[possiblevalues2+1]
  draws <- couplingdeduplication:::coupled_multinomial_(log_N_weights1, log_N_weights2, runif(1), runif(1))
  state1$N <- draws[1] + 1
  state2$N <- draws[2] + 1
  return(list(state1 = state1, state2 = state2))
}

