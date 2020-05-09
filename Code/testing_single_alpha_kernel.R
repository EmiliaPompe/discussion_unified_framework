source('single_kernel_alpha.R')

#--------- testing auxiliary functions
p <- 3
alpha_matrix <- matrix(runif(9), ncol =p, byrow = T)
dim(alpha_matrix)
beta_0 <- rnorm(3, 0)

beta_diff <- alpha_to_beta(alpha_matrix, beta_0)

alpha_2 <- beta_to_alpha(beta_diff, beta_0)
all.equal(alpha_2, alpha_matrix)



#---------- code copied from clearer_single_gibbs
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
V <- V[1:50,1:2]
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
## to compute log-likelihood of the different clusters in the partition
sourceCpp("compute_loglikelihood_all_clusters_all_fields.cpp")
sourceCpp("compute_loglikelihood_all_clusters_one_field.cpp")
sourceCpp("compute_loglikelihood_one_cluster_one_field.cpp")
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
  ## update of N
  ## random walk proposal
  N_rw_stepsize <- 20
  # if (verbose) print("update N")
  Nproposal <- sample(x = (N-N_rw_stepsize):(N+N_rw_stepsize), size = 1)
  if (Nproposal >= partition$ksize){
    target_proposal <- lfactorial(Nproposal) - lfactorial(Nproposal - partition$ksize) - (n+g) * log(Nproposal)
    target_current <- lfactorial(N) - lfactorial(N - partition$ksize) - (n+g) * log(N)
    logu <- log(runif(1))
    if (logu < (target_proposal - target_current)){
      N <- Nproposal
      N_accept <- N_accept + 1
    }
  }
  N
  
# ----------- testing the new function  
  single_beta_0_update(beta_diff,beta_0,  mu_0 = m0, s_0_sq = s0^2, s_sq = s^2)
  l <- 2
  icluster <- 3
 
  single_beta_diff_update(beta_diff_j_l = beta_diff[icluster, l],
                          beta_0_l = beta0[l],
                          l = l,
                          icluster = icluster, 
                          clustering = partition,
                          theta_l = theta[[l]], 
                          V = V,
                          s_sq = s^2,
                          partition_ll = partition_ll,
                          proposal_sd = 0.1)
  
  single_full_alpha_update(beta_diff = beta_diff, alpha = beta_to_alpha(beta_diff, beta0),
                           beta_0 = beta_0, clustering = partition,
                           theta_list = theta, 
                           V = V,
                           partition_ll = partition_ll,
                           mu_0 = m0, 
                           s_0_sq = s0^2,
                           s_sq = s^2,
                           proposal_sd = 0.1,
                           p = p,
                           n = n)
  
  
  
  
  