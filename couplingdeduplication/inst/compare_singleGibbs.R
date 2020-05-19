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
# load("~/discussion_unified_framework/single_gibbs1188f66a9d07a.RData")
# load("~/discussion_unified_framework/authorcode11b164f86a377.RData")


# load("~/discussion_unified_framework/single_gibbs11dcb77b07686.RData")
load("~/discussion_unified_framework/single_gibbs1387d234112f7.RData")

load("~/discussion_unified_framework/authorcode11dea15c19bc0.RData")


single_gibbs_runs[[1]]$theta_history[[1]][1:10,1:10]
authors_runs[[1]]$theta[1:10,1:10]
authors_runs[[1]]$PAR[1:10,1:10]
  
matplot(lapply(single_gibbs_runs, function(x) x$ksize_history[100:5000]) %>% bind_cols(), type = 'l')
matplot(lapply(authors_runs, function(x) x$NZ1[100:5000]) %>% bind_cols(), type = 'l')

matplot(lapply(single_gibbs_runs, function(x) x$N_history[100:5000]) %>% bind_cols(), type = 'l')
matplot(lapply(authors_runs, function(x) x$Npop[100:5000]) %>% bind_cols(), type = 'l')

### comparison of number of clusters in the partition
authors_ksize <- unlist(lapply(authors_runs, function(x) x$NZ1[100:5000]))
table(authors_ksize) / length(authors_ksize)
# 
own_ksize <- unlist(lapply(single_gibbs_runs, function(x) x$ksize_history[100:5000]))
table(own_ksize) / length(own_ksize)
# 
hist(authors_ksize, prob = T)
hist(own_ksize, prob = T, add = T, col = rgb(1,0,0,0.5))

## comparison of the number N
authors_N <- unlist(lapply(authors_runs, function(x) x$Npop[100:5000]))
own_N <- unlist(lapply(single_gibbs_runs, function(x) x$N_history[100:5000]))
hist(authors_N, prob = T, nclass = 30)
hist(own_N, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 30)

authors_theta11 <- unlist(lapply(authors_runs, function(x) x$theta[100:5000,1]))
plot(authors_theta11, type = 'l')

own_theta11 <- unlist(lapply(single_gibbs_runs, function(x) x$theta_history[[1]][100:5000,1]))
plot(own_theta11, type = 'l')

hist(authors_theta11, prob = T, nclass = 30)
hist(own_theta11, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 40)

authors_theta12 <- unlist(lapply(authors_runs, function(x) x$theta[100:5000,2]))
plot(authors_theta12, type = 'l')

own_theta12 <- unlist(lapply(single_gibbs_runs, function(x) x$theta_history[[1]][100:5000,2]))
plot(own_theta12, type = 'l')

hist(authors_theta12, prob = T, nclass = 30)
hist(own_theta12, prob = T, add = T, col = rgb(1,0,0,0.5), nclass = 30)



acf(authors_runs[[1]]$theta[100:5000,1])
acf(single_gibbs_runs[[1]]$theta_history[[1]][100:5000,1])
