V <- V - 1
## initialize parmeters 
### test config
stepsize <- 1000
nMCMC <- 600
nBurn <- nMCMC / 2
nRepeats <- 100
L <- 3
concentration <- 1000
### production config
# nMCMC <- 10000
# nBurn <- 5000
# nRepeats <- 1000


## smaller data config
## work with smaller data
V <- V[1:20,1:4]
## numbers of possibilities for each column
dimV <- dimV[1:4]


dime=as.integer(dimV)
## data size 
n=as.integer(nrow(V))
## number of fields
Hdim=as.integer(ncol(V))
H <- Hdim[1]
## lambda is also eta - labels
lambda <- sample(n, size = n, replace = T) 

cumdime <- 1 +  c(0,  cumsum(dime)) 
## theta - parameter for true record value
## frequencies of categories 
p=as.double(unlist(ALPHA)) ## this is a flat vector, alternatively can store p in an array, but this is what the C file is like
p <- p[1: sum(dime)]
## for small dataset need to change the zeros in p 
# p[p == 0] <- 0.005
# for (l in 1:H){
#   p[1 + cumdime[l] + c(0 : (dime[l] - 1))] <- p[cumdime[l] + c(0 : (dime[l] - 1))]  / sum(p[cumdime[l] + c(0 : (dime[l] - 1))])
# }

## alpha - distorion parameter
a <- matrix(0.01, nrow = n, ncol = H)
## hyperparameter for prior on N 
g <- 1.02
## initial population size 
N1 <- 2500


# find another set of parameters for the coupled gibbs 
lambda1 <- lambda 
lambda2 <- sample(n, size = n, replace = T)
N1 <- N1
N2 <- N1 + 10 
a1 <- matrix(0.01, nrow = n, ncol = H)
a2 <- matrix(0.01, nrow = n, ncol = H)
p1 <- p 
p2 <- rep(NA, length(p1))
for (l in 1:H){
  p2[cumdime[l] : (cumdime[l+1] -1 ) ] <- gtools::rdirichlet(1, alpha = p1[cumdime[l] : (cumdime[l+1] - 1)] * concentration )[1,]
  # print(sum(p2[cumdime[l] : (cumdime[l+1] -1 ) ]))
}

## used to get the ratio
# pcluster_M <- prob_M <- matrix(NA, nrow = n, ncol = H)
# ## before mcmc
# clsize <- rep(0, n)
# for( l in 1:n){
#   clsize[lambda[l] ] <- clsize[lambda[l] ] + 1
# }
# ksize <- 0 
# ksize <- sum(clsize > 0) ## clsize == table(lambda) and ksize = length(clsize)
# psamp <- double(n);
# N <- N1 


### copy and pasted from .c file, will need these later
# am=as.double(rep(log(0.01/0.99),ncol(V)))
# sima=double(nMCMC*ncol(V)) # do not need yet 
# simlambda=integer(nMCMC*nrow(V)) # lambda  
# simnz=integer(nMCMC) # cluster size 
# simNpop=integer(nMCMC) # do not need
# prior=c(0)
# prior=as.integer(prior)
# up.lambda=c(0.15)
# up.lambda=as.double(up.lambda)
# sigma=c(0.5, 0.1, 0.01)
# gzeta=as.double(g)
# sigmal=as.double(sigma)
# simp=double(nMCMC*sum(dimV))
# N1=250
