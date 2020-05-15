## This script aims at coupling the Gibbs sampler on the synthetic data set
## this is work in progress
rm(list = ls())
set.seed(1)
## load packages
library(clusteval)
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
## to perform full sweep of eta updates 
# sourceCpp("update_eta.cpp")
## to perform full sweep of theta updates
source("update_theta.R")
source('single_kernel_alpha.R')
concentration <- 10000 
## define logit and expit transformation
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))

## to perform coupled update of eta
sourceCpp("compute_include_exclude_clustering_likelihood.cpp")
sourceCpp("update_clustering.cpp")
source("couple_multinomial_alt.R")
source("couple_eta.R")
## to perform coupled update of theta
source("couple_dirichlet.R")
source("couple_theta.R")
## to perform coupled update of N
## it uses truncation of N to N_max (define a s a hyper parameter below)
source("couple_N.R")
source('coupled_kernel_alpha.R')
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
## precompute terms so that we only have to look them up
lfactorials_precomp <- sapply(0:N_max, function(NN) lfactorial(NN))
lns_precomp <-  sapply(0:N_max, function(NN) log(NN))

## number of MCMC iterations
nmcmc <- 2e2
## there should be update frequencies ... 

## whether to print some things during the run, or not
verbose <- TRUE

## record history of certain components
N_history1 <- N_history2 <- rep(NA, nmcmc)
# theta1_field1_history <- theta2_field1_history <- matrix(NA, nrow = nmcmc, ncol = Mvec[1])
beta0_history1 <- beta0_history2 <- matrix(NA, nrow = nmcmc, ncol = p)
ksize_history1 <- ksize_history2 <- rep(NA, nmcmc) ## number of clusters

## initialize two chains from different points
## initial etas
eta1 <- sample(n, size = n, replace = T)
# eta2 <- eta1
eta2 <- sample(n, size = n, replace = T)
partition1 <- init_clustering_cpp(eta1 - 1)
partition2 <- init_clustering_cpp(eta2 - 1)
relabel_result1 <- relabel2(eta1, partition1)
relabel_result2 <- relabel2(eta2, partition2)
eta1 <- relabel_result1$eta
partition1$clsize <- partition1$clsize[relabel_result1$old_to_new] 
partition1$clmembers <- partition1$clmembers[relabel_result1$old_to_new,]
eta2 <- relabel_result2$eta
partition2$clsize <- partition2$clsize[relabel_result2$old_to_new] 
partition2$clmembers <- partition2$clmembers[relabel_result2$old_to_new,]
## initial N
N1 <- max(n, 2 * length(unique(eta1)))
N2 <- N1
# N2 <- max(n, 2 * length(unique(eta2)))
# if (N1 == N2) N2 <- N1 + 10

## initial b0
beta0_chain1 <- rnorm(p, m0, s0)
beta0_chain2 <- beta0_chain1
# beta0_chain2 <- rnorm(p, m0, s0)
## initial beta
beta_diff1 <- beta_diff2 <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p){
	beta_diff1[,field] <- rnorm(n, 0, s)
	beta_diff2[,field] <- beta_diff1[,field]
	# beta_diff2[,field] <- rnorm(n, 0, s)
}
beta_diff_identical <- matrix(FALSE, nrow = n, ncol = p)
## compute alphas
alpha1 <- beta_to_alpha(beta_diff1, beta0_chain1)
alpha2 <- beta_to_alpha(beta_diff2, beta0_chain2)
## initial frequencies of categories, initial values
theta1 <- fieldfrequencies[1:p]
## initialize theta2 with a dirichlet move from theta1
theta2 <- fieldfrequencies[1:p]

logtheta1 <- lapply(theta1, function(l) log(l))
logtheta2 <- lapply(theta2, function(l) log(l))

# theta2 <- list()
# for (field in 1:p){
# 	theta2[[field]] <- gtools::rdirichlet(1, alpha = (1 + concentration * theta1[[field]]))[1,]
# }
## sanity check
all(Mvec[1:p] == sapply(theta1, length))
all(Mvec[1:p] == sapply(theta2, length))

## compute log-likelihood associated with each cluster for both chains
partition_ll1 <- compute_loglikelihood_all_clusters_all_fields_cpp(partition1, theta1, logtheta1, V-1, alpha1)
partition_ll2 <- compute_loglikelihood_all_clusters_all_fields_cpp(partition2, theta2, logtheta2, V-1, alpha2)
uponetajoining_loglikelihood1 <- partition_ll1
uponetajoining_loglikelihood2 <- partition_ll2


## initialize quantities to monitor 
cluster_similarity_history <- rep(NA, nmcmc)

## precomputation
s0_sq = s0^2
s_sq = s^2
proposal_sd = sqrt(0.5)

## start the coupled gibbs updates
## note: in the lagged version, we can run some imcmcations of single chain 
## updates on chain1
for (imcmc in 1:nmcmc){
	## update of eta 
  # init_clustering_cpp(eta1-1)$clsize
  # print(partition1$clsize)  
  coupled_update_eta_result <- coupled_update_eta_cpp(eta1-1, eta2-1, partition1, partition2, 
                                                      partition_ll1, partition_ll2, 
                                                      uponetajoining_loglikelihood1, 
                                                      uponetajoining_loglikelihood2, 
                                                      theta1, theta2, logtheta1, logtheta2, 
                                                      V-1, alpha1, alpha2, N1, N2)
  # update_eta_result <- update_eta_cpp(eta1-1, partition1, partition_ll1, uponetajoining_loglikelihood1, theta1, logtheta1, V-1, alpha1, 100000)
  # update_eta_result <- update_eta_cpp(eta2-1, partition2, partition_ll2, uponetajoining_loglikelihood2, theta2, logtheta2, V-1, alpha2, N2)
  
  eta1 <- coupled_update_eta_result$eta1 + 1
  eta2 <- coupled_update_eta_result$eta2 + 1
  partition_ll1 <- coupled_update_eta_result$clusterloglikelihoods1
  partition_ll2 <- coupled_update_eta_result$clusterloglikelihoods2
  partition1 <- coupled_update_eta_result$clustering1
  partition2 <- coupled_update_eta_result$clustering2
  ksize_history1[imcmc] <- partition1$ksize
  ksize_history2[imcmc] <- partition2$ksize
  # print(partition1$clsize)  
  
	## relabel lambd1 and lambda2 and change a1, a2, partiton_ll1, partition_ll2 accordingly
	relabel_result1 <- relabel2(eta1, partition1)
	eta1 <- relabel_result1$eta
	partition1$clsize <- partition1$clsize[relabel_result1$old_to_new] 
	partition1$clmembers <- partition1$clmembers[relabel_result1$old_to_new,]
	partition_ll1 <- partition_ll1[relabel_result1$old_to_new,]
	cl_size1 <- partition1$clsize
	alpha1 <- alpha1[relabel_result1$old_to_new,]
	beta_diff1 <- beta_diff1[relabel_result1$old_to_new,]
	# print(partition1$clsize)  
	#
	relabel_result2 <- relabel2(eta2, partition2)
	eta2 <- relabel_result2$eta 
	partition2$clsize <- partition2$clsize[relabel_result2$old_to_new] 
	partition2$clmembers <- partition2$clmembers[relabel_result2$old_to_new,]
	partition_ll2 <- partition_ll2[relabel_result2$old_to_new,]
	cl_size2 <- partition2$clsize
	alpha2 <- alpha2[relabel_result2$old_to_new,]
	beta_diff2 <- beta_diff2[relabel_result2$old_to_new,]
	#
	ksize_history1[imcmc] <- partition1$ksize 
	ksize_history2[imcmc] <- partition2$ksize 
	cluster_similarity_history[imcmc] <- clusteval::cluster_similarity(eta1, eta2)
	if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'cluster similarity = ', clusteval::cluster_similarity(eta1,eta2), "\n")
	
	### update of N1 and N2
	## using maximal coupling
	log_N_weights1 <- rep(-Inf, N_max)
	possiblevalues1 <- partition1$ksize:N_max
	log_N_weights1[possiblevalues1] <- lfactorials_precomp[possiblevalues1+1] -
	  lfactorials_precomp[possiblevalues1+1-partition1$ksize] -
	  (n+g) * lns_precomp[possiblevalues1+1]
	max_log_N_weights1 <- max(log_N_weights1)
	N_weights1 <- exp(log_N_weights1 - max_log_N_weights1)
	log_N_weights2 <- rep(-Inf, N_max)
	possiblevalues2 <- partition2$ksize:N_max
	log_N_weights2[possiblevalues2] <- lfactorials_precomp[possiblevalues2+1] -
	  lfactorials_precomp[possiblevalues2+1-partition2$ksize] -
	  (n+g) * lns_precomp[possiblevalues2+1]
	max_log_N_weights2 <- max(log_N_weights2)
	N_weights2 <- exp(log_N_weights2 - max_log_N_weights2)
	N1N2 <- couple_multinomial_alt(N_weights1/sum(N_weights1), N_weights2/sum(N_weights2), N_max)
	N1 <- N1N2[1]
	N2 <- N1N2[2]
	N_history1[imcmc] <- N1
	N_history2[imcmc] <- N2

	## update alpha 
	# couple_alpha_result <- coupled_full_alpha_update(beta_diff_1 = beta_diff1,
	#                                                  beta_diff_2 = beta_diff2,
	#                                                  beta_diff_identical = beta_diff_identical,
	#                                                  alpha_1 = alpha1,
	#                                                  alpha_2 = alpha2,
	#                                                  beta_0_1 = beta0_chain1,
	#                                                  beta_0_2 = beta0_chain2,
	#                                                  clustering_1 = partition1,
	#                                                  clustering_2 = partition2,
	#                                                  theta_list_1 = theta1, 
	#                                                  theta_list_2 = theta2,
	#                                                  V = V,
	#                                                  partition_ll_1 = partition_ll1,
	#                                                  partition_ll_2 = partition_ll2,
	#                                                  mu_0 = m0, 
	#                                                  s_0_sq = s0^2,
	#                                                  s_sq = s^2,
	#                                                  proposal_sd = sqrt(0.5),
	#                                                  p = p,
	#                                                  n = n,
	#                                                  cl_size_1 = partition1$clsize,
	#                                                  cl_size_2 = partition2$clsize) 
	# beta0_chain1 <- couple_alpha_result$beta_0_1
	# beta0_chain2 <- couple_alpha_result$beta_0_2
	# beta_0_identical <- couple_alpha_result$beta_0_identical
	# beta_diff1 <- couple_alpha_result$beta_diff_1
	# beta_diff2 <- couple_alpha_result$beta_diff_2
	# beta_diff_identical <- couple_alpha_result$beta_diff_identical
	# alpha1 <- couple_alpha_result$alpha_1
	# alpha2 <- couple_alpha_result$alpha_2
	
	## update theta
	# couple_theta_result <- couple_theta(theta1, theta2, partition1, partition2, partition_ll1, partition_ll2, alpha1, alpha2)
	# theta1 <- couple_theta_result$theta1
	# theta2 <- couple_theta_result$theta2
	# theta1_field1_history[imcmc,] <- theta1[[1]]
	# theta2_field1_history[imcmc,] <- theta2[[1]]
	if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize =', partition1$ksize, '/', partition2$ksize, 'N = ', N1, '/',N2,  "\n")
}

# par(mfrow = c(1,2))
## look at trace and histogram of number of clusters
tail_length = 50
## note: xlim and ylim here works for the full dataset
plot(tail(ksize_history1, tail_length), type = 'l')
lines(tail(ksize_history2, tail_length), type = 'l', col = 'red')

# ksize1_post <- table(tail(ksize1_history, tail_length)) / tail_length
# ksize2_post <- table(tail(ksize2_history, tail_length)) / tail_length
# plot(as.numeric(names(ksize1_post)),ksize1_post,xlim=c(400,500),xlab="k",type="b",lwd=1,cex.lab=2,
#      main="ksize",ylab="", col = 'red')
# points(as.numeric(names(ksize2_post)),ksize2_post,xlim=c(410,495),xlab="k",type="b",lwd=1,cex.lab=2)
# 
# ## look at trace and histogram of population size
plot(tail(N_history1, tail_length), type = 'l')
lines(tail(N_history2, tail_length), type = 'l', col = 'red')
# 
# hist(tail(N1_history, tail_length), xlim = c(500,4500), col=rgb(1,0,0,0.5))
# hist(tail(N2_history, tail_length), col=rgb(0,0,1,0.5), add=T)
