## This script aims at simply running the Gibbs sampler on a small data set
## this is work in progress
rm(list = ls())
set.seed(1)
library(couplingdeduplication)
## load packages
library(Rcpp)
library(ggplot2)
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
hyper$g <- 2
# g <- 1.02

concentration <- 10000
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
lfactorials <- sapply(0:N_max, function(NN) lfactorial(NN))
lns <-  sapply(0:N_max, function(NN) log(NN))
precomp <- list()
precomp$N_max <- N_max
precomp$lfactorials <- lfactorials
precomp$lns <- lns
precomp$proposal_sd = sqrt(0.5)



## whether to print some things during the run, or not
verbose <- TRUE
## number of MCMC iterations
nmcmc <- 1e2
##
single_gibbs_run <- single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, precomp = precomp, verbose = verbose)

## trace plots
par(mfrow = c(1,3))
matplot(single_gibbs_run$ksize_history[1:nmcmc], type = 'l')
matplot(single_gibbs_run$N_history[1:nmcmc], type = "l")
matplot(single_gibbs_run$beta_0_history[1:nmcmc,], type = 'l')
par(mfrow = c(2,2))
for (field in 1:p){
  matplot(single_gibbs_run$theta_history[[field]][1:nmcmc,], type = 'l')
}

# save(N_history, beta_0_history, theta1_history, ksize_history, nmcmc, file = "~/test_single_gibbs.RData")

