## This script aims at simply running the Gibbs sampler on a small data set
## this is work in progress
rm(list = ls())
set.seed(1)
## load packages
library(Rcpp)
library(ggplot2)
## set graphical themes
theme_set(theme_minimal())
## load data 'V' and defines n and p, the dimensions of V
library(RecordLinkage)
data(RLdata500)
myRLDATA=cbind(id=identity.RLdata500,RLdata500)
M1=matrix(unlist(strsplit(soundex(myRLDATA$fname_c1),split="")),ncol=4,byrow=TRUE)
M2=matrix(unlist(strsplit(soundex(myRLDATA$lname_c1),split="")),ncol=4,byrow=TRUE)
M3=matrix(unlist(strsplit(as.character(myRLDATA[,6]),split="")),ncol=4,byrow=TRUE)
M4=matrix( character(nrow(myRLDATA)*2),ncol=2)
for (i in 1:2) M4[,i]=as.character(myRLDATA[,6+i])

V=cbind(M1,M2,M3,M4)
V=as.data.frame(V)
p=ncol(V);

fieldfrequencies=list()
for (l in 1:p) fieldfrequencies[[l]]=table(V[,l])/nrow(V)
Mvec <- c()
for (l in 1:p) {
  Mvec[l]=length(levels(V[,l]))
  V[,l]=as.numeric(V[,l])
}
V=as.matrix(V)
# V is a matrix of n * p
# entries are values in 1:M_l, where M_l is the number of possible categories in field l
# and field l in in 1:p
# Vector 'Mvec' has entries 'M_l' for l in 1:p

## let's subset the data 
V <- V[301:500,1:5]
Mvec <- Mvec[1:5]
fieldfrequencies <- fieldfrequencies[1:5]
## define dimensions of V
n <- dim(V)[1]
p <- dim(V)[2]
g <- 1.1

## get some 'ground truth'
# MATCH=outer(myRLDATA$id, myRLDATA$id, FUN="==")
# MATCH=MATCH[upper.tri(MATCH)]

## source functions
## to relabel 
source("relabel.R")
## to define clustering/partition given vector of labels 'eta'
sourceCpp("init_clustering.cpp")
sourceCpp('compute_loglikelihood.cpp')
## to compute log-likelihood of the different clusters in the partition
#sourceCpp("compute_loglikelihood_all_clusters_all_fields.cpp")
#sourceCpp("compute_loglikelihood_all_clusters_one_field.cpp")
#sourceCpp("compute_loglikelihood_one_cluster_one_field.cpp")
## to perform full sweep of eta updates 
# sourceCpp("update_eta.cpp")
## to perform full sweep of theta updates
source("update_theta.R")
source('single_kernel_alpha.R')
concentration <- 10000 
## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))


## hyper parameter specification
## prior parameter for N 
# g <- 1.02
## prior parameter for beta_0
m0 <- logit(0.01)
s0 <- sqrt(0.1)
## prior parameter for beta given beta_0
s  <- sqrt(0.5)
## truncation limit for N
N_max <- 100000;
## precompute terms so that we only have to look them up
lfactorials_precomp <- sapply(0:N_max, function(NN) lfactorial(NN))
lns_precomp <-  sapply(0:N_max, function(NN) log(NN))



## number of MCMC iterations
nmcmc <- 5e3
## there should be update frequencies ... 

## whether to print some things during the run, or not
verbose <- TRUE

## record history of certain components
N_history <- rep(NA, nmcmc)
theta1_history <- matrix(NA, nrow = nmcmc, ncol = Mvec[1])
beta_0_history <- matrix(NA, nrow = nmcmc, ncol = p)
ksize_history <- rep(NA, nmcmc)

## The state of the Markov chain  is (eta, N, b', b0, theta)

## Initialization of the chains
## initial eta 
eta <- sample(x = 1:n, size = n, replace = T)
## initial N
N <- max(n, 2 * length(unique(eta)))
## 
partition <- init_clustering_cpp(eta-1)
relabel_result <- relabel2(eta, partition)
# clusteval::cluster_similarity(eta, relabel_result$eta)
eta <- relabel_result$eta
partition$clsize <- partition$clsize[relabel_result$old_to_new] 
partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
## initial b0
beta_0 <- rnorm(p, m0, s0)
## initial beta
beta_diff <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p) beta_diff[,field] <- rnorm(n, 0, s)
# beta_diff[partition$clsize==0,] <- NA

##
alpha <- beta_to_alpha(beta_diff, beta_0)
## initial frequencies of categories, initial values
theta <- fieldfrequencies[1:p]
logtheta <- lapply(theta, function(l) log(l))
## sanity check
all(Mvec[1:p] == sapply(theta, length))
## compute log-likelihood associated with each cluster
partition_ll <- compute_loglikelihood_all_clusters_all_fields_cpp(partition, theta, logtheta, V-1, alpha)
uponetajoining_loglikelihood <- partition_ll
## -2.244316 -2.603690 -3.649659
## initialize quantities to monitor
N_accept <- 0
theta_accept <- rep(0, p)

## precomputation
s0_sq = s0^2
s_sq = s^2
proposal_sd = sqrt(0.5)

## next iterate updates in the Gibbs sampler
for (imcmc in 1:nmcmc){
  if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize = ', partition$ksize, 'N = ', N,  "\n")
  ## update of eta
  # print(N)
  # print(length(unique(eta)))
  update_eta_result <- update_eta_cpp(eta-1, partition, partition_ll, uponetajoining_loglikelihood, theta, logtheta, V-1, alpha, N)
  eta <- update_eta_result$eta + 1
  partition_ll <- update_eta_result$clusterloglikelihoods
  partition <- update_eta_result$clustering
  ksize_history[imcmc] <- partition$ksize
  ## relabel 
  relabel_result <- relabel2(eta, partition)
  eta <- relabel_result$eta
  partition$clsize <- partition$clsize[relabel_result$old_to_new]
  partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
  partition_ll <- partition_ll[relabel_result$old_to_new,]
  cl_size <- partition$clsize
  alpha <- alpha[relabel_result$old_to_new,]
  beta_diff <- beta_diff[relabel_result$old_to_new,]
  ## update of N
  ## truncate N to N_max
  log_N_weights <- rep(-Inf, N_max)
  possiblevalues <- partition$ksize:N_max
  log_N_weights[possiblevalues] <- lfactorials_precomp[possiblevalues+1] -
    lfactorials_precomp[possiblevalues+1-partition$ksize] -
    (n+g) * lns_precomp[possiblevalues+1]
  max_log_N_weights <- max(log_N_weights)
  N_weights <- exp(log_N_weights - max_log_N_weights)
  N <- sample.int(n = N_max, size = 1, prob = N_weights)
  N_history[imcmc] <- N
  ## update of  beta_0
  ## using the formula for the conjugate posterior in the Normal-Normal model
  common_var <- 1/(1/s0_sq + n/s_sq)
  for (field in 1:p){
    # calculating beta_prime as previous_beta_0 + new differences
    beta_prime <- beta_0[field] + beta_diff[,field]
    # posterior mean in the Normal-Normal model
    mu <- common_var*(m0/s0_sq + sum(beta_prime)/s_sq)
    beta_0[field] <- rnorm(1, mean = mu, sd = sqrt(common_var))
  }
  # ## alternatively, partially collapsing by dropping the empty clusters
  # common_var <- 1/(1/s0_sq + partition$ksize/s_sq)
  # for (field in 1:p){
  #   # calculating beta_prime as previous_beta_0 + new differences
  #   beta_prime <- beta_0[field] + beta_diff[partition$clsize>0,field]
  #   # posterior mean in the Normal-Normal model
  #   mu <- common_var*(m0/s0_sq + sum(beta_prime)/s_sq)
  #   beta_0[field] <- rnorm(1, mean = mu, sd = sqrt(common_var))
  # }
  beta_0_history[imcmc,] <- beta_0
  ## update of alpha/beta
  for (icluster in 1:n){
    ## for each cluster...
    if (cl_size[icluster] <= 1){
      ## this cluster is empty or contains only one element so likelihood
      ## does not depend on beta parameter; thus beta follows the prior
      beta_diff_cl <- rnorm(p, mean = 0, sd = s)
      exp_beta_ <- exp(beta_0+beta_diff_cl)
      alpha_cl <- exp_beta_/(exp_beta_+1)
      # result_list <- list(beta_diff_j_l = beta_diff_j_l, alpha_j_l = alpha_j_l, accept = NA)
      beta_diff[icluster,] <- beta_diff_cl
      alpha[icluster,] <- alpha_cl
    } else {
      ## cluster has at least 2 elements, we need the likelihood into account
      # loop over fields
      for (l in 1:p){
        ## current beta diff
        current_beta_diff <- beta_diff[icluster, l]
        ## we need to compute the cluster log-likelihood
        proposed_beta_diff <- rnorm(1, current_beta_diff, sd = proposal_sd/max(1,cl_size[icluster]))
        ## calculating the corresponding value of alpha
        exp_beta <- exp(beta_0[l] + proposed_beta_diff)
        proposed_alpha <- exp_beta/(exp_beta+1)
        ## random number for acceptance/rejection
        logu <- log(runif(1))
        ## log likelihood taken from partition_ll + the prior
        current_posterior <-  partition_ll[icluster, l] + dnorm(current_beta_diff, 0, s, log = TRUE)
        ## we use beta_diff hence the mean is 0, not beta_0
        proposed_logprior <- dnorm(proposed_beta_diff, 0, s, log = TRUE)
        proposed_loglik <- compute_loglikelihood_one_cluster_one_field_cpp(l-1, icluster-1, partition,
                                                                           theta[[l]], logtheta[[l]], V-1, proposed_alpha)
        accept <- (logu < (proposed_loglik + proposed_logprior - current_posterior))
        if (accept){
          ## update log-likelihood associated cluster
          partition_ll[icluster,l] <- proposed_loglik
          beta_diff[icluster,l] <- proposed_beta_diff
          alpha[icluster,l] <- proposed_alpha
        } else {
          ## other wise do nothing
        }
      }
    }
  }
  ## update of theta
  ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
  update_theta_result <- update_theta(theta, partition, partition_ll, alpha)
  theta <- update_theta_result$theta
  logtheta <- lapply(theta, function(l) log(l))
  theta_accept <- theta_accept + update_theta_result$theta_accept
  theta1_history[imcmc, ] <- theta[[1]]
  partition_ll <- update_theta_result$partition_ll
}

# cat(N_accept/nmcmc, "\n")
# cat(theta_accept/nmcmc, "\n")

# par(mfrow = c(2,1))
nburn <- floor(nmcmc / 2)
## reproduce top-left plot of figure 2
matplot(ksize_history[nburn : nmcmc], type = 'l')
plot(N_history, type = "l")

matplot(beta_0_history, type = 'l')

# matplot(theta1_history, type = 'l')

hist(beta_0_history[5e2:nmcmc,1])
# save(beta_0_history, file = "~/beta0noncollapsed.RData")
# save(beta_0_history, file = "~/beta0partiallycollapsed.RData")

# save(N_history, beta_0_history, theta1_history, ksize_history, nmcmc, file = "~/cleanergibbs_run1.RData")

## prior of ksize
# rzeta=function(a){
#   b=2^{a-1}
#   cond=TRUE
#   while(cond){
#     u=runif(1)
#     v=runif(1)
#     x=floor(u^(-1/(a-1)))
#     if (x==Inf) x=.Machine$double.xmax
#     t=(1+1/x)^(a-1)
#     cond=(v*x*(t-1)/(b-1)>t/b)
#   }
#   return(x)
# }
# NN2=c()
# for(i in 1:100000) NN2[i]=rzeta(1.02)
# k2=c()
# for(i in 1:100000) {
#   if (NN2[i]>.Machine$integer) k2[i]=length(unique(sample(.Machine$integer,size=500,rep=T)))
#   else k2[i]=length(unique(sample(NN2[i],size=500,rep=T)))
# }
# ksize_prior=table(k2)/length(k2)
# ksize_prior_names=as.numeric(names(ksize_prior))
# lines(ksize_prior_names,ksize_prior,col=2,lwd=2)
# 
# ## reproduce the top-right plot of figure 2
# matplot(N_history[nburn : nmcmc], type = 'l')
# points(N_history[nburn : nmcmc])
# hist(N_history[nburn : nmcmc])
# 
# 
# 
# matplot(theta1_history[,1:2], type = 'l')
# abline(h = fieldfrequencies[[1]][1:2])
# 
# # pairs(theta1_history[200:nmcmc,1:10])
# # hist(N_history[(nmcmc/10):nmcmc])
# 
# hist(N_history[(nmcmc/10):nmcmc])
# hist(theta1_history[(nmcmc/10):nmcmc,1])
# abline(v = fieldfrequencies[[1]][1])
# 
# hist(theta1_history[(nmcmc/10):nmcmc,2])
# abline(v = fieldfrequencies[[1]][2])
# 
# hist(theta1_history[(nmcmc/10):nmcmc,3])
# abline(v = fieldfrequencies[[1]][3])
# 
# # density plots for beta_0
# ggplot(as.data.frame(beta_0_history), aes(V1)) + geom_density()
# ggplot(as.data.frame(beta_0_history), aes(V2)) + geom_density()
# ggplot(as.data.frame(beta_0_history), aes(V14)) + geom_density()
# 
# # traceplots for beta_0
# matplot(beta_0_history[seq(from = 1, by = 1, to = 1e4),14], type = 'l')
# 
# p1 <- ggplot(data.frame(V1= N_history[2001:nmcmc]), aes(V1)) + geom_histogram(bins=20, col = 'black', fill = 'gray70') + xlab('N')
# p2 <- ggplot(data.frame(V1= N_history[1:nmcmc], V2 = 1:nmcmc), aes(x=V2,y =V1)) + geom_line(col = 'black')+ ylab('N') + xlab('iteration')
# p1
# p2
# library(gridExtra)
# grid.arrange(p1, p2, ncol =2)
# #ggsave('single_chain_N_hist_traceplot.png', grid.arrange(p1, p2, ncol =2), height = 6, width =10)

