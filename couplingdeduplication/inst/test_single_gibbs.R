## This script aims at simply running the Gibbs sampler on a small data set
## this is work in progress
rm(list = ls())
set.seed(1)
library(couplingdeduplication)
## load packages
library(ggplot2)
library(doParallel)
library(doRNG)
library(dplyr)
registerDoParallel(cores = detectCores()-2)
## set graphical themes
theme_set(theme_minimal())

RLdata <- get_RLdata500()
## let's subset the data 
# selectedrows <- 401:420
# selectedcolumns <- 4:7
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
# hyper$g <- 2
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
precomp <- list()
precomp$N_max <- N_max
precomp$lfactorials <- lfactorials
precomp$lns <- lns
precomp$proposal_sd = sqrt(0.5)
precomp$concentration <- 10000





# single_gibbs_run <- single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, precomp = precomp, verbose = verbose)
single_gibbs <- function(nmcmc, V, fieldfrequencies, hyper, precomp, verbose = TRUE){
  n <- dim(V)[1]
  p <- dim(V)[2]
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## draw initial state
  state <- rinit(n, fieldfrequencies, hyper)
  state$alpha <- matrix(0.01, nrow = n, ncol = p)
  state$beta_0 <- NA
  state$beta_diff <- matrix(log(0.01/0.99), nrow = n, ncol = p)
  ## record history of certain components
  N_history <- rep(NA, nmcmc)
  ksize_history <- rep(NA, nmcmc)
  theta_history <- list()
  for (field in 1:p){
    theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field])
  }
  
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
    ## update of theta
    ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
    update_theta_result <- update_theta(state$theta, state$partition, partition_ll, state$alpha, precomp$concentration)
    state$theta <- update_theta_result$theta
    state$logtheta <- lapply(state$theta, function(x) log(x))
    for (field in 1:p){
      theta_history[[field]][imcmc, ] <- state$theta[[field]]
    }
    partition_ll <- update_theta_result$partition_ll
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
  }
  return(list(ksize_history = ksize_history, 
              N_history = N_history, theta_history = theta_history))
}

## whether to print some things during the run, or not
verbose <- TRUE
## number of MCMC iterations
nmcmc <- 5e3
##
# single_gibbs_run <- single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, precomp = precomp, verbose = verbose)
# plot(single_gibbs_run$ksize_history, type = 'l')
# plot(single_gibbs_run$N_history, type = 'l')
# matplot(single_gibbs_run$theta_history[[1]][,1:3], type = 'l')

single_gibbs_runs <- foreach(irep = 1:6) %dorng% {
  single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, precomp = precomp, verbose = verbose)
}

# par(mfrow = c(1,1))
# matplot(lapply(single_gibbs_runs, function(x) x$ksize_history) %>% bind_cols(), type = 'l')
# matplot(lapply(single_gibbs_runs, function(x) x$N_history) %>% bind_cols(), type = 'l')

filename_ <- tempfile(pattern = "single_gibbs", tmpdir = "~/discussion_unified_framework", fileext = ".RData")
save(single_gibbs_runs, nmcmc, n, p, V, fieldfrequencies, hyper, precomp, file = filename_)

matplot(lapply(single_gibbs_runs, function(x) x$ksize_history) %>% bind_cols(), type = 'l')
matplot(lapply(single_gibbs_runs, function(x) x$N_history) %>% bind_cols(), type = 'l')


# ## trace plots
# par(mfrow = c(1,3))
# single_gibbs_run <- single_gibbs_runs[[1]]
# matplot(single_gibbs_run$ksize_history[1:nmcmc], type = 'l')
# matplot(single_gibbs_run$N_history[1:nmcmc], type = "l")
# matplot(single_gibbs_run$beta_0_history[1:nmcmc,], type = 'l')
# par(mfrow = c(2,2))
# for (field in 1:p){
#   matplot(single_gibbs_run$theta_history[[field]][1:nmcmc,], type = 'l')
# }

