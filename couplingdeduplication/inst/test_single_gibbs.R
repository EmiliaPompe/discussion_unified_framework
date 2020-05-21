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
algotuning$eta_update_prob <- 0.1
algotuning$verbose <- TRUE
## number of MCMC iterations
nmcmc <- 1e2
## single run
single_gibbs_run <- single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies,
                                 hyper = hyper, algotuning = algotuning, update.theta = TRUE)
##
plot(single_gibbs_run$ksize_history[10:nmcmc], type = 'l')
plot(single_gibbs_run$N_history[10:nmcmc], type = 'l')
matplot(single_gibbs_run$theta_history[[1]][,1:3], type = 'l')



## multiple runs, in parallel
# library(doParallel)
# library(doRNG)
# registerDoParallel(cores = detectCores()-2)
# single_gibbs_runs <- foreach(irep = 1:6) %dorng% {
#   single_gibbs(nmcmc = nmcmc, V = V, fieldfrequencies = fieldfrequencies, hyper = hyper, algotuning = algotuning,
#                update.theta = FALSE, verbose = verbose)
# }

## save output in RData file
## filename_ <- tempfile(pattern = "single_gibbs_update_eta_N_theta", tmpdir = "~/discussion_unified_framework", fileext = ".RData")
## filename_ <- "~/discussion_unified_framework/single_gibbs_update_eta_N.RData"
## save(single_gibbs_runs, nmcmc, n, p, V, fieldfrequencies, hyper, algotuning, file = filename_)
## save(single_gibbs_runs, nmcmc, n, p, V, fieldfrequencies, hyper, algotuning, file = filename_)
