#'@title coupled_gibbs
#'@description  runs two chains, with a lag of at least one, and until either meeting occurs, or iteration 'm',
#' and at most 'max_iterations'
#' at the moment this only couples the update of eta and N
#'@param V a matrix representing the data
#'@param fieldfrequencies a list ...
#'@param hyper list containing hyperparameters
#'@param algotuning list containing precomputed and tuning parameters
#'@param m time horizon: the chains are sampled until the maximum between m and the meeting time; equal to one
#'@param lag time lag, equal to one by default
#'@param max_iterations maximum number of iterations, at which to interrup the while loop; Inf by default
#'@export
coupled_gibbs <- function(V, fieldfrequencies, hyper, algotuning, update.theta = FALSE, m = 1, lag = 1, max_iterations = Inf){
  starttime <- Sys.time()
  ## dimension of the data
  ## number of rows
  n <- dim(V)[1]
  ## number of fields
  p <- dim(V)[2]
  ## number of possibilities in each field
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## current time
  time <- 0
  ## record history of certain components
  preallocate <- 1e2
  nrowsamples1 <- m+preallocate+lag
  N_history1 <- rep(NA, nrowsamples1)
  ksize_history1 <- rep(NA, nrowsamples1)
  theta_history1 <- list()
  N_history2 <- rep(NA, nrowsamples1-lag)
  ksize_history2 <- rep(NA, nrowsamples1-lag)
  theta_history2 <- list()
  for (field in 1:p){ 
    theta_history1[[field]] <- matrix(NA, nrow = nrowsamples1, ncol = Mvec[field]) 
    theta_history2[[field]] <- matrix(NA, nrow = nrowsamples1-lag, ncol = Mvec[field]) 
  }
  ## draw initial state number one 
  state1 <- rinit(n, fieldfrequencies, hyper, V)
  N_history1[time+1] <- state1$N
  ksize_history1[time+1] <- state1$partition$ksize
  for (field in 1:p){
    theta_history1[[field]][time+1,] <- state1$theta[[field]]
  }
  ## advance first chain by 'lag' iterations
  if (algotuning$verbose){  cat("advancing first chain for", lag, "iterations...\n") }
  for (iteration in 1:lag){
    time <- time + 1
    ## update eta and relabel
    state1 <- update_eta_relabel(state1, V, algotuning)
    ## update N given rest
    state1 <- update_N(state1, V, hyper, algotuning)
    ##
    if (update.theta){
      state1 <- update_theta(state1, V, algotuning)
    }
    N_history1[time+1] <- state1$N
    ksize_history1[time+1] <- state1$partition$ksize
    for (field in 1:p){
      theta_history1[[field]][time+1,] <- state1$theta[[field]]
    }
  }
  if (algotuning$verbose){  cat("...done. Now, coupled updates... \n") }
  ## draw second initial state
  state2 <- rinit(n, fieldfrequencies, hyper, V)
  N_history2[1] <- state2$N
  ksize_history2[1] <- state2$partition$ksize
  for (field in 1:p){
    theta_history2[[field]][1,] <- state2$theta[[field]]
  }
  ## meeting time
  meetingtime <- Inf
  ## next iterate updates in the coupled Gibbs sampler
  while ((time < max(meetingtime, m)) && (time < max_iterations)){
    time <- time + 1 # time is lag+1,lag+2,...
    if (algotuning$verbose && (time %% 10 == 0)){
      cat("iteration", time, '\n chain1: ksize = ', state1$partition$ksize, 'N = ', state1$N,
          '\n chain2: ksize = ', state2$partition$ksize, 'N = ', state2$N,  "\n")
      cat("cluster similarity between chains", clusteval::cluster_similarity(state1$eta, state2$eta), "\n")
    }
    if (is.finite(meetingtime)){
      state1 <- update_eta_relabel(state1, V, algotuning)
      state1 <- update_N(state1, V, hyper, algotuning)
      if (update.theta){
        state1 <- update_theta(state1, V, algotuning)
      }
      state2 <- state1
    } else {
      coupled_update_eta_relabel_results <- coupled_update_eta_relabel(state1, state2, V, algotuning)
      state1 <- coupled_update_eta_relabel_results$state1
      state2 <- coupled_update_eta_relabel_results$state2
      ## indicator that all eta variables are equal
      equalcomponents <- all(state1$eta == state2$eta)
      ## update of N
      ## truncate N to N_max
      results_ <- coupled_update_N(state1, state2, V, hyper, algotuning)
      state1 <- results_$state1
      state2 <- results_$state2
      equalcomponents <- equalcomponents && all(state1$N == state2$N)
      ## update of theta
      if (update.theta){
        results_ <- coupled_update_theta(state1, state2, V, algotuning)
        state1 <- results_$state1
        state2 <- results_$state2
        equalcomponents <- equalcomponents && results_$identical
      }
      if (is.infinite(meetingtime) && equalcomponents){
        meetingtime <- time
        if (algotuning$verbose){
          cat("!! meeting at time", time, "\n")
          cat("distance between theta:\n")
          print(sapply(1:p, function(field) sum(abs(state1$theta[[field]] - state2$theta[[field]]))))
          cat("distance between N:\n")
          print(state1$N - state2$N)
          cat("distance between ksize:\n")
          print(state1$partition$ksize - state2$partition$ksize)
          cat("partition ll:\n")
          print(sapply(1:p, function(field) sum(abs(state1$partition_ll[,field] - state2$partition_ll[,field]), na.rm = T)))
        }
      }
    }
    if ((time+1) > nrowsamples1){
      new_rows <- nrowsamples1
      nrowsamples1 <- nrowsamples1 + new_rows
      ksize_history1 <- c(ksize_history1, rep(NA, new_rows))
      ksize_history2 <- c(ksize_history2, rep(NA, new_rows))
      N_history1 <- c(N_history1, rep(NA, new_rows))
      N_history2 <- c(N_history2, rep(NA, new_rows))
      for (field in 1:p){
        theta_history1[[field]] <- rbind(theta_history1[[field]], matrix(NA, nrow = new_rows, ncol = Mvec[field]))
        theta_history2[[field]] <- rbind(theta_history2[[field]], matrix(NA, nrow = new_rows, ncol = Mvec[field]))
      }
    }
    N_history1[time+1] <- state1$N
    ksize_history1[time+1] <- state1$partition$ksize
    N_history2[time+1-lag] <- state2$N
    ksize_history2[time+1-lag] <- state2$partition$ksize
    for (field in 1:p){
      theta_history1[[field]][time+1,]     <- state1$theta[[field]]
      theta_history2[[field]][time+1-lag,] <- state2$theta[[field]]
    }
  }
  ksize_history1 <- ksize_history1[1:(time+1)]
  ksize_history2 <- ksize_history2[1:(time-lag+1)]
  N_history1 <- N_history1[1:(time+1)]
  N_history2 <- N_history2[1:(time-lag+1)]
  for (field in 1:p){
    theta_history1[[field]] <- theta_history1[[field]][1:(time+1),]
    theta_history2[[field]] <- theta_history2[[field]][1:(time-lag+1),]
  }
  cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  if (algotuning$verbose){ cat("elapsed time for this run:", elapsedtime, "seconds \n") }
  return(list(ksize_history1 = ksize_history1, 
              ksize_history2 = ksize_history2, 
              N_history1 = N_history1,
              N_history2 = N_history2, 
              theta_history1 = theta_history1, 
              theta_history2 = theta_history2,
              meetingtime = meetingtime,
              cost = cost, elapsedtime = elapsedtime))
}
