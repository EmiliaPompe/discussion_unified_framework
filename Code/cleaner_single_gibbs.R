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
#V <- V[1:50,1:2]
## define dimensions of V
n <- dim(V)[1]
p <- dim(V)[2]

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
sourceCpp("update_eta.cpp")
## to perform full sweep of theta updates
source("update_theta.R")
source('single_kernel_alpha.R')
concentration <- 10000 
## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))


## hyper parameter specification
## prior parameter for N 
g <- 1.02
## prior parameter for beta0
m0 <- logit(0.01)
s0 <- sqrt(0.1)
## prior parameter for beta given beta0
s  <- sqrt(0.5)
## truncation limit for N
N_max <- 100000;

## number of MCMC iterations
nmcmc <- 1e2
## there should be update frequencies ... 

## whether to print some things during the run, or not
verbose <- TRUE

## record history of certain components
N_history <- rep(NA, nmcmc)
theta1_history <- matrix(NA, nrow = nmcmc, ncol = Mvec[1])
beta0_history <- matrix(NA, nrow = nmcmc, ncol = p)

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
beta0 <- rnorm(p, m0, s0)
## initial beta
beta_diff <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p) beta_diff[,field] <- rnorm(n, 0, s)
# beta_diff[partition$clsize==0,] <- NA

##
alpha <- beta_to_alpha(beta_diff, beta0)
## initial frequencies of categories, initial values
theta <- fieldfrequencies[1:p]
## sanity check
all(Mvec[1:p] == sapply(theta, length))
## compute log-likelihood associated with each cluster
partition_ll <- compute_loglikelihood_all_clusters_all_fields_cpp(partition, theta, V-1, alpha)
## initialize quantities to monitor
N_accept <- 0
theta_accept <- rep(0, p)
## next iterate updates in the Gibbs sampler
for (imcmc in 1:nmcmc){
  if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, "\n")
  ## update of eta
  # print(N)
  # print(length(unique(eta)))
  update_eta_result <- update_eta_cpp(eta-1, partition, partition_ll, theta, V-1, Mvec[1:p], alpha, N, beta0, s)
  eta <- update_eta_result$eta
  partition_ll <- update_eta_result$clusterloglikelihoods
  partition <- update_eta_result$clustering
  # alpha <- update_eta_result$alpha
  ## relabel 
  relabel_result <- relabel2(eta, partition)
  eta <- relabel_result$eta
  partition$clsize <- partition$clsize[relabel_result$old_to_new]
  partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
  partition_ll <- partition_ll[relabel_result$old_to_new,]
  alpha <- alpha[relabel_result$old_to_new,]
  beta_diff <- beta_diff[relabel_result$old_to_new,] 
  ## update of N

  # ## random walk proposal
  # N_rw_stepsize <- 20
  # # if (verbose) print("update N")
  # Nproposal <- sample(x = (N-N_rw_stepsize):(N+N_rw_stepsize), size = 1)
  # if (Nproposal >= partition$ksize){
  #   target_proposal <- lfactorial(Nproposal) - lfactorial(Nproposal - partition$ksize) - (n+g) * log(Nproposal)
  #   target_current <- lfactorial(N) - lfactorial(N - partition$ksize) - (n+g) * log(N)
  #   logu <- log(runif(1))
  #   if (logu < (target_proposal - target_current)){
  #     N <- Nproposal
  #     N_accept <- N_accept + 1
  #   }
  # }
  # N_history[imcmc] <- N

  ## truncate N to N_max
  # if (verbose) print("update N")
  log_N_weights <- sapply(1 : N_max, function(NN) lfactorial(NN) - lfactorial(NN - partition$ksize) - (n+g) * log(NN) )
  max_log_N_weights <- max(log_N_weights)
  log_N_weights <- log_N_weights - max_log_N_weights
  N_weights <- exp(log_N_weights) / sum(exp(log_N_weights))
  N <- sample.int(n = N_max, size = 1, prob = N_weights)
  N_history[imcmc] <- N
  # if (verbose) print("update N done")
  
  # update of alpha, beta diff and beta0
  # we are storing both alpha and beta diff to avoid making unnecessary calculations
  alpha_beta_update <- single_full_alpha_update(beta_diff = beta_diff, 
                                                alpha = alpha,
                                                beta_0 = beta0, 
                                                clustering = partition,
                                                theta_list = theta, 
                                                V = V,
                                                partition_ll = partition_ll,
                                                mu_0 = m0, 
                                                s_0_sq = s0^2,
                                                s_sq = s^2,
                                                proposal_sd = 1.5,
                                                p = p,
                                                n = n)
  alpha <- alpha_beta_update$alpha
  beta_diff <- alpha_beta_update$beta_diff
  beta0 <- alpha_beta_update$beta0
  beta0_history[imcmc,] <- beta0
  print(mean(alpha_beta_update$acc_matrix, na.rm = TRUE))
 
  
  ## update of theta
  ## note: prior on theta = uniform on simplex, equivalently Dirichlet(1,1,...,1)
  update_theta_result <- update_theta(theta, partition, partition_ll, alpha)
  theta <- update_theta_result$theta
  theta_accept <- theta_accept + update_theta_result$theta_accept
  theta1_history[imcmc, ] <- theta[[1]]
  partition_ll <- update_theta_result$partition_ll
}
 
cat(N_accept/nmcmc, "\n")
cat(theta_accept/nmcmc, "\n")

matplot(N_history, type = 'l')

matplot(theta1_history[,1:2], type = 'l')
abline(h = fieldfrequencies[[1]][1:2])

# pairs(theta1_history[200:nmcmc,1:10])
# hist(N_history[(nmcmc/10):nmcmc])

hist(N_history[(nmcmc/10):nmcmc])
hist(theta1_history[(nmcmc/10):nmcmc,1])
abline(v = fieldfrequencies[[1]][1])

hist(theta1_history[(nmcmc/10):nmcmc,2])
abline(v = fieldfrequencies[[1]][2])

hist(theta1_history[(nmcmc/10):nmcmc,3])
abline(v = fieldfrequencies[[1]][3])

# density plots for beta0
ggplot(as.data.frame(beta0_history), aes(V1)) + geom_density()
ggplot(as.data.frame(beta0_history), aes(V2)) + geom_density()
ggplot(as.data.frame(beta0_history), aes(V14)) + geom_density()

# traceplots for beta0
matplot(beta0_history[seq(from = 1, by = 1, to = 1e4),14], type = 'l')

p1 <- ggplot(data.frame(V1= N_history[2001:nmcmc]), aes(V1)) + geom_histogram(bins=20, col = 'black', fill = 'gray70') + xlab('N')
p2 <- ggplot(data.frame(V1= N_history[1:nmcmc], V2 = 1:nmcmc), aes(x=V2,y =V1)) + geom_line(col = 'black')+ ylab('N') + xlab('iteration')
p1
p2
library(gridExtra)
grid.arrange(p1, p2, ncol =2)
ggsave('single_chain_N_hist_traceplot.png', grid.arrange(p1, p2, ncol =2), height = 6, width =10)

