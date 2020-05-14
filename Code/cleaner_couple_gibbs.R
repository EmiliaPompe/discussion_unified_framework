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
# V <- V[1:50,1:2]
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
## to perform full sweep of eta updates 
sourceCpp("update_eta.cpp")
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

## number of MCMC imcmcations
nmcmc <- 1e2
## there should be update frequencies ... 

## whether to print some things during the run, or not
verbose <- TRUE

## record history of certain components
N1_history <- N2_history <- rep(NA, nmcmc)
theta1_field1_history <- theta2_field1_history <- matrix(NA, nrow = nmcmc, ncol = Mvec[1])
beta0_chain1_history <- beta0_chain2_history <- matrix(NA, nrow = nmcmc, ncol = p)
ksize1_history <- ksize2_history <- rep(NA, nmcmc) ## number of clusters

## initialize two chains from different points
## initial etas
eta1 <- sample(n, size = n, replace = T)
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
N2 <- max(n, 2 * length(unique(eta2)))
if (N1 == N2) N2 <- N1 + 10
## initial b0
beta0_chain1 <- rnorm(p, m0, s0)
beta0_chain2 <- rnorm(p, m0, s0)
## initial beta
beta_diff1 <- beta_diff2 <- matrix(NA, nrow = n, ncol = p)
for (field in 1:p){
	beta_diff1[,field] <- rnorm(n, 0, s)
	beta_diff2[,field] <- rnorm(n, 0, s)
}
## compute alphas
alpha1 <- beta_to_alpha(beta_diff1, beta0_chain1)
alpha2 <- beta_to_alpha(beta_diff2, beta0_chain2)
## initial frequencies of categories, initial values
theta1 <- fieldfrequencies[1:p]
## initialize theta2 with a dirichlet move from theta1
theta2 <- list()
for (field in 1:p){
	theta2[[field]] <- gtools::rdirichlet(1, alpha = (1 + concentration * theta1[[field]]))[1,]
}
## sanity check
all(Mvec[1:p] == sapply(theta1, length))
all(Mvec[1:p] == sapply(theta2, length))

## compute log-likelihood associated with each cluster for both chains
partition_ll1 <- compute_loglikelihood_all_clusters_all_fields_cpp(partition1, theta1, V-1, alpha1)
partition_ll2 <- compute_loglikelihood_all_clusters_all_fields_cpp(partition2, theta1, V-1, alpha2)

## intailize quantities to monitor 
cluster_distance_history <- rep(NA, nmcmc)
## start the coupled gibbs updates
## note: in the lagged version, we can run some imcmcations of single chain 
## updates on chain1
for (imcmc in 1 : nmcmc){
	## update of eta 
	for (j in 1:n){
	  couple_eta_result <- couple_eta(j, eta1, eta2, partition1, partition2, partition_ll1, partition_ll2, theta1, theta2, alpha1, alpha2, N1, N2)
	  eta1 <- couple_eta_result$eta1
	  eta2 <- couple_eta_result $eta2
	  partition1 <- couple_eta_result$partition1
	  partition2 <- couple_eta_result$partition2
	  partition_ll1 <- couple_eta_result$partition_ll1
	  partition_ll2 <- couple_eta_result$partition_ll2
	}
	## relabel lambd1 and lambda2 and change a1, a2, partiton_ll1, partition_ll2 accordingly
	relabel_result1 <- relabel2(eta1, partition1)
	eta1 <- relabel_result1$eta
	partition1$clsize <- partition1$clsize[relabel_result1$old_to_new] 
	partition1$clmembers <- partition1$clmembers[relabel_result1$old_to_new,]
	partition_ll1 <- partition_ll1[relabel_result1$old_to_new,]
	alpha1 <- alpha1[relabel_result1$old_to_new,]
	beta_diff1 <- beta_diff1[relabel_result1$old_to_new,]
	relabel_result2 <- relabel2(eta2, partition2)
	eta2 <- relabel_result2$eta
	partition2$clsize <- partition2$clsize[relabel_result2$old_to_new] 
	partition2$clmembers <- partition2$clmembers[relabel_result2$old_to_new,]
	partition_ll2 <- partition_ll2[relabel_result2$old_to_new, ]
	alpha2 <- alpha2[relabel_result2$old_to_new,]
	beta_diff2 <- beta_diff2[relabel_result2$old_to_new,]
	ksize1_history[imcmc] <- partition1$ksize 
	ksize2_history[imcmc] <- partition2$ksize 
	cluster_distance_history[imcmc] <- clusteval::cluster_similarity(eta1,eta2)
	if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'cluster distance = ', clusteval::cluster_similarity(eta1,eta2), "\n")
	## update alpha 
	## to be added 
	
	## update theta
	couple_theta_result <- couple_theta(theta1, theta2, partition1, partition2, partition_ll1, partition_ll2, alpha1, alpha2)
	theta1 <- couple_theta_result$theta1
	theta2 <- couple_theta_result$theta2
	theta1_field1_history[imcmc,] <- theta1[[1]]
	theta2_field1_history[imcmc,] <- theta2[[1]]
	## update of N 
	Ns <- couple_N(N1 = N1, N2 = N2, k1 = partition1$ksize, k2 = partition2$ksize, g, n)
	N1 <- Ns[1]
	N2 <- Ns[2]
	N1_history[imcmc] <- N1
	N2_history[imcmc] <- N2
	if (verbose && (imcmc %% 1 == 0)) cat("iteration", imcmc, "/", nmcmc, 'ksize =', partition1$ksize, '/', partition2$ksize, 'N = ', N1, '/',N2,  "\n")
}

par(mfrow = c(1,2))
## look at trace and histogram of number of clusters
tail_length = 50
## note: xlim and ylim here works for the full dataset
plot(tail(ksize1_history, tail_length), type = 'l', ylim = c(400,500))
lines(tail(ksize1_history, tail_length), type = 'l', col = 'red') 
ksize1_post <- table(tail(ksize1_history, tail_length)) / tail_length
ksize2_post <- table(tail(ksize2_history, tail_length)) / tail_length
plot(as.numeric(names(ksize1_post)),ksize1_post,xlim=c(400,500),xlab="k",type="b",lwd=1,cex.lab=2,
     main="ksize",ylab="", col = 'red')
points(as.numeric(names(ksize2_post)),ksize2_post,xlim=c(410,495),xlab="k",type="b",lwd=1,cex.lab=2)

## look at trace and histogram of population size
plot(tail(N1_history, tail_length), type = 'l', ylim = c(500,4500))
lines(tail(N2_history, tail_length), type = 'l', col = 'red')

hist(tail(N1_history, tail_length), xlim = c(500,4500), col=rgb(1,0,0,0.5))
hist(tail(N2_history, tail_length), col=rgb(0,0,1,0.5), add=T)
