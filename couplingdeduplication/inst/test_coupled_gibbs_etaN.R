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
precomp <- list()
precomp$N_max <- N_max
precomp$lfactorials <- lfactorials
precomp$lns <- lns
precomp$proposal_sd = sqrt(0.5)
precomp$concentration <- 10000

nmcmc <- 1e2

coupled_gibbs <- function(nmcmc, V, fieldfrequencies, hyper, precomp, verbose = TRUE){
  n <- dim(V)[1]
  p <- dim(V)[2]
  Mvec <- unlist(lapply(fieldfrequencies, function(l) length(l)))
  ## draw initial state
  state1 <- rinit(n, fieldfrequencies, hyper)
  state1$alpha <- matrix(0.01, nrow = n, ncol = p)
  state1$beta_0 <- NA
  state1$beta_diff <- matrix(log(0.01/0.99), nrow = n, ncol = p)
  state2 <- rinit(n, fieldfrequencies, hyper)
  state2$alpha <- matrix(0.01, nrow = n, ncol = p)
  state2$beta_0 <- NA
  state2$beta_diff <- matrix(log(0.01/0.99), nrow = n, ncol = p)
  
  ## record history of certain components
  N_history1 <- rep(NA, nmcmc)
  ksize_history1 <- rep(NA, nmcmc)
  N_history2 <- rep(NA, nmcmc)
  ksize_history2 <- rep(NA, nmcmc)
  
  ## compute log-likelihood associated with each cluster
  partition_ll1 <- couplingdeduplication:::compute_loglikelihood_all_clusters_all_fields_cpp(state1$partition, state1$theta, 
                                                                                             state1$logtheta, V-1, state1$alpha)
  partition_ll2 <- couplingdeduplication:::compute_loglikelihood_all_clusters_all_fields_cpp(state2$partition, state2$theta, 
                                                                                             state2$logtheta, V-1, state2$alpha)
  
  meetingtime <- Inf
  
  ## next iterate updates in the coupled Gibbs sampler
  for (imcmc in 1:nmcmc){
    if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize = ', state1$partition$ksize, 'N = ', state1$N,  "\n")
    
    update_eta_result <- couplingdeduplication:::coupled_update_eta(state1$eta-1, state2$eta-1, state1$partition, state2$partition, 
                                                                    partition_ll1, partition_ll2, state1$theta, state2$theta, 
                                                                    state1$logtheta, state2$logtheta, V-1, 
                                                                    state1$alpha, state2$alpha, state1$N, state2$N, 0.1)
    state1$eta <- update_eta_result$eta1 + 1
    partition_ll1 <- update_eta_result$clusterloglikelihoods1
    state1$partition <- update_eta_result$clustering1
    ksize_history1[imcmc] <- state1$partition$ksize
    state2$eta <- update_eta_result$eta2 + 1
    partition_ll2 <- update_eta_result$clusterloglikelihoods2
    state2$partition <- update_eta_result$clustering2
    ksize_history2[imcmc] <- state2$partition$ksize
    ## 
    if (verbose) cat("NA in partition ll make sense?", all(is.na(partition_ll1[,1]) == (state1$partition$clsize==0)), 
                     all(is.na(partition_ll2[,1]) == (state2$partition$clsize==0)), "\n")
    
    ## relabel 
    relabel_result1 <- relabel(state1$eta, state1$partition)
    relabel_result2 <- relabel(state2$eta, state2$partition)
    #
    state1$eta <- relabel_result1$eta
    state1$partition$clsize <- state1$partition$clsize[relabel_result1$old_to_new]
    state1$partition$clmembers <- state1$partition$clmembers[relabel_result1$old_to_new,]
    partition_ll1 <- partition_ll1[relabel_result1$old_to_new,]
    state1$alpha <- state1$alpha[relabel_result1$old_to_new,]
    state1$beta_diff <- state1$beta_diff[relabel_result1$old_to_new,]
    
    state2$eta <- relabel_result2$eta
    state2$partition$clsize <- state2$partition$clsize[relabel_result2$old_to_new]
    state2$partition$clmembers <- state2$partition$clmembers[relabel_result2$old_to_new,]
    partition_ll2 <- partition_ll2[relabel_result2$old_to_new,]
    state2$alpha <- state2$alpha[relabel_result2$old_to_new,]
    state2$beta_diff <- state2$beta_diff[relabel_result2$old_to_new,]
    ## indicator that all eta variables are equal
    eta_equal <- all(state1$eta == state2$eta)
    #
    # clusteval::cluster_similarity(state1$eta, state2$eta)
    # partition_ll1 - partition_ll2 
    # 
    # ## update of theta
    # ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
    # ## for each field
    # for (field in 1:p){
    #   concentration <- precomp$concentration
    #   ## current state
    #   x <- state$theta[[field]] 
    #   cl_log_lik <- partition_ll[,field] 
    #   ## dirichlet proposal
    #   x_propose <- rgamma(n = length(x), shape = concentration*x+1, rate = 1)
    #   x_propose <- x_propose/sum(x_propose)
    #   logx_propose <- log(x_propose)
    #   cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
    #   ## transition ratio 
    #   lratio <- sum((concentration*x_propose) * log(x) - (concentration*x) * log(x_propose)) + sum(lgamma(concentration*x+1) - lgamma(concentration*x_propose+1))
    #   ## log accept probability
    #   laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
    #   if (log(runif(1)) < laccept){
    #     state$theta[[field]] <- x_propose
    #     state$logtheta[[field]] <- log(x_propose)
    #     partition_ll[,field] <- cl_log_lik_new
    #   }
    #   theta_history[[field]][imcmc, ] <- state$theta[[field]]
    # }
    ## update of N
    ## truncate N to N_max
    log_N_weights1 <- rep(-Inf, precomp$N_max)
    log_N_weights2 <- rep(-Inf, precomp$N_max)
    possiblevalues1 <- (state1$partition$ksize):(precomp$N_max)
    possiblevalues2 <- (state2$partition$ksize):(precomp$N_max)
    log_N_weights1[possiblevalues1] <- precomp$lfactorials[possiblevalues1+1] -
      precomp$lfactorials[possiblevalues1+1-state1$partition$ksize] -
      (state1$n+hyper$g) * precomp$lns[possiblevalues1+1]
    log_N_weights2[possiblevalues2] <- precomp$lfactorials[possiblevalues2+1] -
      precomp$lfactorials[possiblevalues2+1-state2$partition$ksize] -
      (state2$n+hyper$g) * precomp$lns[possiblevalues2+1]
    draws <- couplingdeduplication:::coupled_multinomial_(log_N_weights1, log_N_weights2, runif(1), runif(1))
    state1$N <- draws[1] + 1
    state2$N <- draws[2] + 1
    N_history1[imcmc] <- state1$N
    N_history2[imcmc] <- state2$N
    N_equal <- all(state1$N == state2$N)
    if (is.infinite(meetingtime) && N_equal && eta_equal){
      meetingtime <- imcmc
    }
    if (is.finite(meetingtime)){
      if (!N_equal || !eta_equal){
        print("chains unmet")
      }
    }
  }
  return(list(ksize_history1 = ksize_history1, 
              ksize_history2 = ksize_history2, 
              N_history1 = N_history1,
              N_history2 = N_history2, 
              meetingtime = meetingtime))
}

# draws <- t(sapply(1:1e4, function(i) couplingdeduplication:::coupled_multinomial_(log_N_weights1, log_N_weights2, runif(1), runif(1))))
# plot(table(draws[,1]) / 1e4)
# w1 <- exp(log_N_weights1 - max(log_N_weights1))
# w1 <- w1 / sum(w1)
# hist(sample(x = seq_along(w1), size = 1e4, prob = w1, replace = T), prob = T, add = T, col = 'red', nclass = 100)
# plot(table(draws[,2]) / 1e4)
# w2 <- exp(log_N_weights2 - max(log_N_weights2))
# w2 <- w2 / sum(w2)
# hist(sample(x = seq_along(w2), size = 1e4, prob = w2, replace = T), prob = T, add = T, col = 'red', nclass = 100)



## whether to print some things during the run, or not
verbose <- TRUE
## number of MCMC iterations
nmcmc <- 1e2
##

coupled_gibbs_run <- coupled_gibbs(nmcmc, V, fieldfrequencies, hyper, precomp, verbose)
coupled_gibbs_run$meetingtime
plot(coupled_gibbs_run$ksize_history1, type = 'l')
lines(coupled_gibbs_run$ksize_history2, col = 'red')
abline(v = coupled_gibbs_run$meetingtime)
plot(coupled_gibbs_run$N_history1, type = 'l')
lines(coupled_gibbs_run$N_history2, type = 'l', col = 'red')
abline(v = coupled_gibbs_run$meetingtime)

