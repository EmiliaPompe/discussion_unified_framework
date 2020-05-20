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

