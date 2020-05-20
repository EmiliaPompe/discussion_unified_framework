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
algotuning$theta_update_conc <- 10000
algotuning$eta_update_prob <- 0.1
algotuning$verbose <- TRUE


lag <- 50
nrep <- 50
### uncomment following code to obtain some meeting times
# coupled_gibbs_runs <- foreach(irep = 1:nrep) %dorng% {
#   algotuning$verbose <- FALSE
#   coupled_gibbs(V, fieldfrequencies, hyper, algotuning, m = 1, lag = lag, max_iterations = 1e4)
# }
# save(lag, nrep, coupled_gibbs_runs, file = "coupledgibbs_update_eta_N.RData")
## load results
load(file = "coupledgibbs_update_eta_N.RData")
## obtain meeting times
meeting_times <- sapply(coupled_gibbs_runs, function(x) x$meetingtime)
## obtain TV upper bounds
tv_upper_bound_estimates <- function(meeting_times, L, t){
  return(mean(pmax(0,ceiling((meeting_times-L-t)/L))))
}
## plot TV upper bounds
niter <- (floor(1.1*max(meeting_times)-lag))
upperbounds <- sapply(1:niter, function(t) tv_upper_bound_estimates(unlist(meeting_times), lag, t))
hist_noplot <- hist(meeting_times - lag, plot = F, nclass = 30)
xgrid <- c(min(hist_noplot$breaks), hist_noplot$mids, max(hist_noplot$breaks))
densitygrid <- c(0, hist_noplot$density, 0)
g_tvbounds <- qplot(x = 1:niter, y = upperbounds, geom = "line")
g_tvbounds <- g_tvbounds + ylab("TV upper bounds") + xlab("iteration")
g_tvbounds <- g_tvbounds + scale_y_continuous(breaks = c(0,1), limits = c(0,1.1))
g_tvbounds <- g_tvbounds + geom_ribbon(data=data.frame(x = xgrid,
                                                       ymin = rep(0, length(xgrid)),
                                                       y = densitygrid/max(densitygrid)),
                                       aes(x= x, ymin = ymin, ymax = y, y=NULL), alpha = .3) + geom_line()
print(g_tvbounds)
