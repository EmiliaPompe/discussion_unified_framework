## This script aims at experimenting with coupled Gibbs samplers 
## this is work in progress
rm(list = ls())
set.seed(1)
library(couplingdeduplication)
## load packages
library(ggplot2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
## set graphical themes
theme_set(theme_minimal())

RLdata <- get_RLdata500()
selectedrows <- 1:nrow(RLdata$V)
selectedcolumns <- 1:ncol(RLdata$V)
V <- RLdata$V[selectedrows, selectedcolumns]
Mvec <- RLdata$Mvec[selectedcolumns]
fieldfrequencies <- RLdata$fieldfrequencies[selectedcolumns]
## define dimensions of V
n <- dim(V)[1]
p <- dim(V)[2]
## hyper parameter specification
hyper <- list()
## prior parameter on N
hyper$g <- 1.02

## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))

## prior parameter for beta_0
hyper$m0 <- logit(0.01)
hyper$s0 <- sqrt(0.1)
## prior parameter for beta given beta_0
hyper$s  <- sqrt(0.5)
## some precomputation
hyper$s0_sq = hyper$s0^2
hyper$s_sq = hyper$s^2
## truncation limit for N
N_max <- 1e5;
## precompute terms so that we only have to look them up
## and store algorithmic tuning 
lfactorials <- sapply(0:N_max, function(NN) lfactorial(NN))
lns <-  sapply(0:N_max, function(NN) log(NN))
algotuning <- list()
algotuning$N_max <- N_max
algotuning$lfactorials <- lfactorials
algotuning$lns <- lns
algotuning$proposal_sd = sqrt(0.5)
algotuning$theta_update_rwconcentration <- 10000
algotuning$theta_update_rw <- 0.5
algotuning$theta_update_indepscale <- 0.8
algotuning$eta_update_prob <- 0.1
algotuning$verbose <- TRUE

## coupled Gibbs without the re-labeling step
## to see the effect on the meeting time
coupled_gibbs <- function(V, fieldfrequencies, hyper, algotuning, m = 1, lag = 1, max_iterations = Inf){
  starttime <- Sys.time()
  ## dimension of the data
  ## number of rows
  n <- dim(V)[1]
  ## number of fields
  p <- dim(V)[2]
  ## number of possibilities in each field
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  time <- 0
  ## record history of certain components
  preallocate <- 1e2
  nrowsamples1 <- m+preallocate+lag
  N_history1 <- rep(NA, nrowsamples1)
  ksize_history1 <- rep(NA, nrowsamples1)
  N_history2 <- rep(NA, nrowsamples1-lag)
  ksize_history2 <- rep(NA, nrowsamples1-lag)
  ## 
  ## draw initial state number one 
  state1 <- rinit(n, fieldfrequencies, hyper, V)
  N_history1[time+1] <- state1$N
  # ksize_history1[time+1] <- state1$partition$ksize
  ## advance first chain by 'lag' iterations
  if (algotuning$verbose){  cat("advancing first chain for", lag, "iterations...\n") }
  for (iteration in 1:lag){
    time <- time + 1
    ## update eta and relabel
    state1 <- update_eta_relabel(state1, V, algotuning)
    ## update N given rest
    state1 <- update_N(state1, V, hyper, algotuning)
    ##
    N_history1[time+1] <- state1$N
    # ksize_history1[time+1] <- state1$partition$ksize
  }
  if (algotuning$verbose){  cat("...done. Now, coupled updates... \n") }
  ## draw second initial state
  state2 <- rinit(n, fieldfrequencies, hyper, V)
  N_history2[1] <- state2$N
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
      state2 <- state1
    } else {
      ## draw new etas from max coupling
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
      ## indicator that all eta variables are equal
      eta_equal <- all(state1$eta == state2$eta)
      ## update of N
      ## truncate N to N_max
      results_ <- coupled_update_N(state1, state2, V, hyper, algotuning)
      state1 <- results_$state1
      state2 <- results_$state2
      N_equal <- all(state1$N == state2$N)
      if (is.infinite(meetingtime) && N_equal && eta_equal){
        meetingtime <- time
        if (algotuning$verbose){ cat("!! meeting at time", time, "\n") }
      }
    }
    if ((time+1) > nrowsamples1){
      new_rows <- nrowsamples1
      nrowsamples1 <- nrowsamples1 + new_rows
      ksize_history1 <- c(ksize_history1, rep(NA, new_rows))
      ksize_history2 <- c(ksize_history2, rep(NA, new_rows))
      N_history1 <- c(N_history1, rep(NA, new_rows))
      N_history2 <- c(N_history2, rep(NA, new_rows))
    }
    
    ksize_history1[time+1] <- state1$partition$ksize
    ksize_history2[time+1-lag] <- state2$partition$ksize
    N_history1[time+1] <- state1$N
    N_history2[time+1-lag] <- state2$N
    
  }
  ksize_history1 <- ksize_history1[1:(time+1)]
  ksize_history2 <- ksize_history2[1:(time-lag+1)]
  N_history1 <- N_history1[1:(time+1)]
  N_history2 <- N_history2[1:(time-lag+1)]
  cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  if (algotuning$verbose){ cat("elapsed time for this run:", elapsedtime, "seconds \n") }
  return(list(ksize_history1 = ksize_history1, 
              ksize_history2 = ksize_history2, 
              N_history1 = N_history1,
              N_history2 = N_history2, 
              meetingtime = meetingtime,
              cost = cost, elapsedtime = elapsedtime))
}




# lag <- 50
# nrep <- 10
# ### uncomment following code to obtain some meeting times
# coupled_gibbs_runs <- foreach(irep = 1:nrep) %dorng% {
#   algotuning$verbose <- FALSE
#   coupled_gibbs(V, fieldfrequencies, hyper, algotuning, m = 1, lag = lag, max_iterations = 1e3)
# }
# save(lag, nrep, coupled_gibbs_runs, file = "coupledgibbs_update_eta_N_norelabel.RData")
## load results
load(file = "coupledgibbs_update_eta_N_norelabel.RData")
names(coupled_gibbs_runs[[1]])
## 
## show traces of one run
coupled_gibbs_run <- coupled_gibbs_runs[[1]]
niter <- length(coupled_gibbs_run$ksize_history1)
plot(0:(niter-1), coupled_gibbs_run$ksize_history1, type = 'l', ylim = c(430, 460))
lines(lag:(niter-1), coupled_gibbs_run$ksize_history2, col = 'red')

plot(0:(niter-1), coupled_gibbs_run$N_history1, type = 'l', ylim = c(1000, 4000))
lines(lag:(niter-1), coupled_gibbs_run$N_history2, col = 'red')

## obtain meeting times
meeting_times <- sapply(coupled_gibbs_runs, function(x) x$meetingtime)
print(meeting_times)

# ## obtain TV upper bounds
# tv_upper_bound_estimates <- function(meeting_times, L, t){
#   return(mean(pmax(0,ceiling((meeting_times-L-t)/L))))
# }
# ## plot TV upper bounds
# niter <- (floor(1.1*max(meeting_times)-lag))
# upperbounds <- sapply(1:niter, function(t) tv_upper_bound_estimates(unlist(meeting_times), lag, t))
# hist_noplot <- hist(meeting_times - lag, plot = F, nclass = 30)
# xgrid <- c(min(hist_noplot$breaks), hist_noplot$mids, max(hist_noplot$breaks))
# densitygrid <- c(0, hist_noplot$density, 0)
# g_tvbounds <- qplot(x = 1:niter, y = upperbounds, geom = "line")
# g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
# g_tvbounds <- g_tvbounds + scale_y_continuous(breaks = c(0,1), limits = c(0,1.1))
# g_tvbounds <- g_tvbounds + geom_ribbon(data=data.frame(x = xgrid,
#                                                        ymin = rep(0, length(xgrid)),
#                                                        y = densitygrid/max(densitygrid)),
#                                        aes(x= x, ymin = ymin, ymax = y, y=NULL), alpha = .3) + geom_line()
# print(g_tvbounds)
