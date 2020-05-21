## This script aims at simply running the Gibbs sampler and seeing if it works
## this is work in progress
rm(list = ls())
set.seed(1)
library(couplingdeduplication)
## load packages
library(ggplot2)
library(dplyr)
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
algotuning$theta_update_conc <- 10000
algotuning$eta_update_prob <- 0.1
algotuning$verbose <- TRUE

## number of MCMC iterations
nmcmc <- 1e2
## single run
# single_gibbs_run <- single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies,
#                                  hyper = hyper, algotuning = algotuning, 
update.theta <- TRUE
##
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
par(mfrow = c(3,1))
plot(ksize_history[10:nmcmc], type = 'l')
plot(N_history[10:nmcmc], type = 'l')
matplot(theta_history[[1]][,1:3], type = 'l')

## now let's do updates of only theta, to see how quickly it mixes when other 
## components are fixed
nmcmc <- 5000
theta_history <- list()
field <- 1
# for (field in 1:p){ 
theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field])
for (imcmc in 1:nmcmc){
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

par(mfrow = c(1,1))
matplot(theta_history[[field]], type = 'l')
acf(theta_history[[field]][,1])
## the mixing is pretty terrible 
## calculate acceptance rate approximately
1-mean(abs(diff(theta_history[[field]][1000:nmcmc,1])<1e-20))

## what if we do MH steps with independent proposals from a Dirichlet distribution?
hist(theta_history[[field]][1000:nmcmc,1], prob = TRUE, nclass = 50)
datacounts <- RLdata$fieldfrequencies[[field]] * state$partition$ksize
post_ <- gtools::rdirichlet(1e5, alpha = datacounts + 1)
hist(post_[,1], add = T, col = rgb(1,0,0,0.5), prob = T, nclass = 50)

nmcmc <- 50000
theta_history <- list()
field <- 1
# for (field in 1:p){ 
theta_history[[field]] <- matrix(NA, nrow = nmcmc, ncol = Mvec[field])
for (imcmc in 1:nmcmc){
  ## with some probability, perform random walk move
  u_ <- runif(1)
  if (u_ < 0.5){
    theta_update_conc <- algotuning$theta_update_conc
    ## current state
    x <- state$theta[[field]] 
    cl_log_lik <- state$partition_ll[,field] 
    ## dirichlet random walk proposal
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
  } else {
    ## else perform independent proposal move 
    alpha_proposal <- RLdata$fieldfrequencies[[field]] * state$partition$ksize * 0.9 + 1
    x_propose <- gtools::rdirichlet(1, alpha = alpha_proposal)
    logx_propose <- log(x_propose)
    cl_log_lik_new <- couplingdeduplication:::compute_loglikelihood_all_clusters_one_field_cpp(field - 1, state$partition, x_propose, logx_propose, V - 1, state$alpha)
    ## current state
    x <- state$theta[[field]] 
    cl_log_lik <- state$partition_ll[,field] 
    ## transition ratio 
    lratio <- log(gtools::ddirichlet(x = x, alpha = alpha_proposal)) - log(gtools::ddirichlet(x = x_propose, alpha = alpha_proposal)) 
    ## log accept probability
    laccept <- lratio + sum(cl_log_lik_new[which(state$partition$clsize != 0)]) - sum(cl_log_lik[which(state$partition$clsize != 0)])
    if (log(runif(1)) < laccept){
      state$theta[[field]] <- x_propose
      state$logtheta[[field]] <- log(x_propose)
      state$partition_ll[,field] <- cl_log_lik_new
    }
    
  }
  theta_history[[field]][imcmc, ] <- state$theta[[field]]
}


par(mfrow = c(1,1))
matplot(theta_history[[field]], type = 'l')
acf(theta_history[[field]][,1])
