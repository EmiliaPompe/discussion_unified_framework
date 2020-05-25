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
algotuning$fieldfrequencies <- fieldfrequencies
algotuning$proposal_sd = sqrt(0.5)
algotuning$theta_update_rwconcentration <- 10000
algotuning$theta_update_rw <- 0.5
algotuning$theta_update_indepscale <- 0.8
algotuning$eta_update_prob <- 0.1
algotuning$verbose <- TRUE
## number of MCMC iterations
nmcmc <- 5e1
## draw initial state
state <- rinit(n, fieldfrequencies, hyper, V)
## next iterate updates in the Gibbs sampler
for (imcmc in 1:nmcmc){
  ## update eta given the rest
  state <- update_eta_relabel(state, V, algotuning)
  ## update of beta
  state <- update_beta(state, hyper, V, algotuning)
  ## update of theta
  state <- update_theta(state, V, algotuning)
  ## update N given rest
  state <- update_N(state, V, hyper, algotuning)
}

## now update only beta

state <- state
state_collapsed <- state
nmcmc <- 20000
beta_0_history <- matrix(NA, nrow = nmcmc, ncol = p)
beta_0_history2 <- matrix(NA, nrow = nmcmc, ncol = p)
for (imcmc in 1:nmcmc){
  state <- update_beta(state, hyper, V, algotuning)
  state_collapsed <- update_beta_collapsed(state_collapsed, hyper, V, algotuning)
  beta_0_history[imcmc,] <- state$beta_0
  beta_0_history2[imcmc,] <- state_collapsed$beta_0
}

matplot(beta_0_history[,1], type = 'l')
matplot(beta_0_history2[,1], type = 'l')

hist(beta_0_history[,1], prob = TRUE, col = 'red', nclass = 100)
hist(beta_0_history2[,1], prob = TRUE, col = rgb(0,0,0,0.5), add= T, nclass = 100)

# 
# acf(beta_0_history[,1])
# acf(beta_0_history2[,1])
