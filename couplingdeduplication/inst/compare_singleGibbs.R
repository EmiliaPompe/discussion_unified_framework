## Script to compare results produced by this package 
## to the results obtained with the authors' code

## Requires having run other scripts before; this script just does some plots
## and calculations to make sure the samplers mostly agree

rm(list = ls())
set.seed(1)
library(couplingdeduplication)
## load packages
library(ggplot2)
## set graphical themes
theme_set(theme_minimal())
setwd("~/discussion_unified_framework/")
## load results from our own sampler

## update only of N and eta

load("~/discussion_unified_framework/single_gibbs_update_eta_N.RData")
# load("~/discussion_unified_framework/single_gibbs_update_eta_N_theta185d1775198c5.RData")
# load("~/discussion_unified_framework/authorcode11b164f86a377.RData")

if ("theta_history" %in% names(single_gibbs_runs[[1]])){
  print(head(single_gibbs_runs[[1]]$theta_history[[1]][,1:5]))
}

sg_nmcmc <- length(single_gibbs_runs[[1]]$ksize_history)
sg_mcmc_subset <- floor(seq(from = sg_nmcmc/10, to = sg_nmcmc, length.out = 1e3))


# load("~/discussion_unified_framework/authorcode11dea15c19bc0.RData")
load("~/discussion_unified_framework/authorcode_lambda_N.RData")

# single_gibbs_runs[[1]]$theta_history[[1]][1:10,1:10]

## is this a run with a theta update?
authors_runs[[1]]$theta[1:10,1:10]
## is this a run with a beta update?
authors_runs[[1]]$PAR[1:10,1:10]


authors_run_nmcmc <- length(authors_runs[[1]]$NZ1)
authors_run_mcmc_subset <- floor(seq(from = authors_run_nmcmc/10, to = authors_run_nmcmc, length.out = 1e3))

## size of the partition, i.e. number of non-empty clusters
matplot(sg_mcmc_subset, lapply(single_gibbs_runs, function(x) x$ksize_history[sg_mcmc_subset]) %>% bind_cols(), type = 'l')
matplot(x = authors_run_mcmc_subset, y = lapply(authors_runs, function(x) x$NZ1[authors_run_mcmc_subset]) %>% bind_cols(), type = 'l')
authors_ksize <- unlist(lapply(authors_runs, function(x) x$NZ1[100:authors_run_nmcmc]))
round(table(authors_ksize) / length(authors_ksize),3)
# 
own_ksize <- unlist(lapply(single_gibbs_runs, function(x) x$ksize_history[100:sg_nmcmc]))
round(table(own_ksize) / length(own_ksize), 3)
hist(authors_ksize, prob = T)
hist(own_ksize, prob = T, add = T, col = rgb(1,0,0,0.5))

## N 
matplot(lapply(single_gibbs_runs, function(x) x$N_history[sg_mcmc_subset]) %>% bind_cols(), type = 'l', ylim = c(1500,4e3))
matplot(lapply(authors_runs, function(x) x$Npop[authors_run_mcmc_subset]) %>% bind_cols(), type = 'l', ylim = c(1500,4e3))

## comparison of the number N
authors_N <- unlist(lapply(authors_runs, function(x) x$Npop[100:authors_run_nmcmc]))
own_N <- unlist(lapply(single_gibbs_runs, function(x) x$N_history[100:sg_nmcmc]))
hist(authors_N, prob = T, nclass = 30)
hist(own_N, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 30)

##
# 
# authors_theta11 <- unlist(lapply(authors_runs, function(x) x$theta[100:5000,1]))
# plot(authors_theta11, type = 'l')
# 
# own_theta11 <- unlist(lapply(single_gibbs_runs, function(x) x$theta_history[[1]][100:10000,1]))
# plot(own_theta11, type = 'l')
# 
# hist(authors_theta11, prob = T, nclass = 30)
# hist(own_theta11, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 50)
# 
# authors_theta12 <- unlist(lapply(authors_runs, function(x) x$theta[100:5000,2]))
# plot(authors_theta12, type = 'l')
# 
# own_theta12 <- unlist(lapply(single_gibbs_runs, function(x) x$theta_history[[1]][100:10000,2]))
# plot(own_theta12, type = 'l')
# 
# hist(authors_theta12, prob = T, nclass = 30)
# hist(own_theta12, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 40)
# 
# acf(authors_runs[[1]]$theta[100:5000,1])
# acf(single_gibbs_runs[[1]]$theta_history[[1]][100:10000,1])
