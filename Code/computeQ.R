rm(list = ls())
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

n <- 20
V <- V[1:n,]
d=ncol(V);


ALPHA=list()
for (l in 1:d) ALPHA[[l]]=table(V[,l])/nrow(V)
dimV=c()
for ( l in 1:d) {
  dimV[l]=length(levels(V[,l]))
  V[,l]=as.numeric(V[,l])
}

V=as.matrix(V)
H <- ncol(V)

MATCH=outer(myRLDATA$id, myRLDATA$id, FUN="==")
MATCH=MATCH[upper.tri(MATCH)]

nMCMC <- 100
## initialization
n <- dim(V)[1]
lambda <- sample(n, size = n, replace = T) 
p=as.double(unlist(ALPHA))
p[p==0] <- 0.005
V <- V - 1
dime=as.integer(dimV)
a=matrix(0.01, nrow = n, ncol = H)
am=as.double(rep(log(0.01/0.99),ncol(V)))
Hdim=as.integer(ncol(V))
nMCMC=as.integer(nMCMC)
sima=double(nMCMC*ncol(V)) # do not need yet 
simlambda=integer(nMCMC*nrow(V)) # lambda  
simnz=integer(nMCMC) # cluster size 
simNpop=integer(nMCMC) # do not need
n=as.integer(nrow(V))
prior=c(0)
prior=as.integer(prior)
up.lambda=c(0.15)
up.lambda=as.double(up.lambda)
g=1.02
sigma=c(0.5, 0.1, 0.01)
gzeta=as.double(g)
sigmal=as.double(sigma)
simp=double(nMCMC*sum(dimV))
N1=2500

## storage
cumdime <- 1 + c(0 ,  cumsum(dime)) ## helps getting the index of theta
H <- Hdim[1]
pcluster_M <- prob_M <- matrix(NA, nrow = n, ncol = H)
## before mcmc
clsize <- rep(0, n)
for( l in 1:n){
  clsize[lambda[l] ] <- clsize[lambda[l] ] + 1
}
ksize <- 0 
ksize <- sum(clsize > 0) ## clsize == table(lambda) and ksize = length(clsize)
psamp <- double(n);
N <- N1

## pcluster_M[z,l] <- marginal of observed cluster
for (j in 1:n){
  if (clsize[lambda[j] ] == 1) ksize <- ksize - 1;
  for (q in 1:n){
    psamp[q] <- 1.0
    old <- F
    old <- any(lambda[-j] == q) ## if old == TRUE, then q labels an observed cluster
    if (old){
      for (h in 1:H){
        # calculate calculate pcluster_M : probability of the cluster q with the xclision of record j  
        pcluster_M[q,h] <- 1.0
        first_element <- TRUE
        for ( i in 1:n){
          if (i != j & lambda[i] == q){
            # print(c(i, first_element))
            if (first_element){
              first_element <- FALSE
              pcluster_M[q,h] <- p[cumdime[h] + V[i,h] ]
            }else{
              cumproduct <- 1.0
              for (k in 1:(i-1)){
                if (k != j & lambda[k] == q) cumproduct <- cumproduct * (( 1- a[q,h]) * (V[k,h] == V[i,h]) + a[q,h] * p[cumdime[h] +  V[k,h]])
              }
              pcluster_M[q,h] <- pcluster_M[q,h] * a[q,h] * p[cumdime[h] + V[i,h]] + (1 - a[q,h]) * p[cumdime[h] + V[i,h]] * cumproduct 
            }
          }
        }
        # calculate numerator 
        prob_M[q, h] <- 1.0
        for (i in 1 : n){
          if (i != j & lambda[i] == q){
            prob_M[q,h] <- prob_M[q,h]*( (V[i,h] == V[j,h])*(1 - a[q,h]) + a[q,h] * p[cumdime[h]+V[i,h]] );
          }
        }
        prob_M[q,h] <- a[q,h]*p[cumdime[h] + V[j,h]] + (1 - a[q,h]) * prob_M[q,h] / pcluster_M[q,h];
        psamp[q] <- psamp[q]*prob_M[q,h];
        if(is.na(psamp[q])) print(c(j,q,h,k,i))
      }
    }else{
      for (h in 1:H){
        psamp[q] = psamp[q] * p[cumdime[h] + V[j,h]]
      }
    }
    if( clsize[q] == 0 | (clsize[q] ==1 & q == lambda[j])) psamp[q] = psamp[q] * (N - ksize ) / (n - ksize)
  }
  ## process psamp 
  psamp <- psamp / sum(psamp)
  clsize[lambda[j]] <- clsize[lambda[j]] - 1 
  lambda[j] <- sample.int(n = n, size = 1, prob = psamp)
  clsize[lambda[j]] <- clsize[lambda[j] ] + 1 
  if (clsize[lambda[j] ] == 1) ksize = ksize + 1
}
