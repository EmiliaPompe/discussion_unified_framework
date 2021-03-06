## This script aims at running this package's implementation of the Gibbs sampler 
## it does so for fixed theta, fixed beta, and then for fixed beta
## (the update of beta is not implemented yet)
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
nmcmc <- 2e2
print("update of eta and N")
gibbs_runs <- foreach(irep = 1:5) %dorng% {
  single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, algotuning = algotuning, update.theta = FALSE)
}
save(gibbs_runs, nmcmc, n, p, V, fieldfrequencies, hyper, algotuning, file = "gibbs_update_eta_N.RData")

print("update of eta, theta and N")
gibbs_runs <- foreach(irep = 1:5) %dorng% {
  single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, algotuning = algotuning, update.theta = TRUE)
}
save(gibbs_runs, nmcmc, n, p, V, fieldfrequencies, hyper, algotuning, file = "gibbs_update_eta_theta_N.RData")

