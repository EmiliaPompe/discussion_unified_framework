#### This file reproduces the two panels of Figure 2 in the paper
#### coupled updates of lambda and N, keeping  theta fixed
#### it does not update alpha, but alpha changes in each iteration due to re-labeling
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 20-April 2020

rm(list = ls())
source("getData.R")
source("initPara.R")
source("coupleMultinomial.R")
library(Rcpp)
sourceCpp("computeNweights.cpp")
sourceCpp("computeQcpp.cpp")
source("coupleLambda.R")
source("coupleN.R")
source("relabel.R")

stepsize <- 1000
nMCMC <- 5000
nRepeats <- 1000
# nMCMC <- 10

coupleGibbs <- function(nMCMC, lambda1, lambda2, N1, N2, g, p, a1, a2, n, earlyStop = TRUE){
  meetTime <- Inf
  N1_chain <- N2_chain <- rep(0, nMCMC)
  k1_chain <- k2_chain <- rep(0, nMCMC)
  percentageCoupled <- rep(0, nMCMC)
  iter <- 1
  for (iter in 1:nMCMC){
    ## update lambda 
    for (j in 1:n){
      lambdajs <- coupleLambda(j, lambda1, lambda2, p1 = p, p2 = p, a1 = a1, a2 = a2, N1 = N1, N2 = N2)
      # if(lambdajs[1] == lambdajs[2]){
      #   print(cat("for j = ", j ," coupled lambda-j = ", lambdajs[1]))
      # }
      lambda1[j] <- lambdajs[1]
      lambda2[j] <- lambdajs[2]
    }
    percentageCoupled[iter] <- mean(lambda1 == lambda2)
    ## relabel lambd1 and lambda2 and change a1 a2 accordingly
    relabel1 <- relabel(lambda1)
    lambda1 <- relabel1$lambda
    a1 <- a1[relabel1$iis,]
    relabel2 <- relabel(lambda2)
    lambda2 <- relabel2$lambda
    a2 <- a2[relabel2$iis,]
    if (iter %% (nMCMC / 50) == 0 ){
      print(paste('iteration', iter, percentageCoupled[iter] * 100, "percent coupled",sep = ' '))
    }
    ksize1 <- length(unique(lambda1))
    ksize2 <- length(unique(lambda2))
    k1_chain[iter] <- ksize1
    k2_chain[iter] <- ksize2
    ## update N 
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize2, g, n)
    N1 <- Ns[1]
    N2 <- Ns[2]
    N1_chain[iter] <- N1
    N2_chain[iter] <- N2
    meet <- FALSE
    if (N1 == N2 & all(lambda1 == lambda2) & !meet){
      meetTime <- iter
      meet <- TRUE
    }
    if (earlyStop & meet) break
  }
  return(list(N1_chain = N1_chain, N2_chain = N2_chain,
              k1_chain = k1_chain, k2_chain = k2_chain,
              meetTime = meetTime, percentageCoupled = percentageCoupled))
}

result <- coupleGibbs(nMCMC, lambda1, lambda2, N1, N2, g, p, a1, a2,n, earlyStop = FALSE)
print(result$meetTime)
print(which(result$N1_chain != result$N2_chain))


# pdf("reFigure2.pdf")
par(mfrow = c(1,2))
## save this plot
hist(result$k1_chain)
hist(result$N1_chain)
# dev.off()
par(mfrow = c(1,2))
plot(result$N1_chain, type = 'l')
lines(result$N2_chain, col = 'red')
plot(result$k1_chain, type = 'l')
lines(result$k2_chain, col = 'red')

pdf("percentageCoupledLambda.pdf")
par(mfrow = c(1,1))
plot(result$percentageCoupled, type  = 'l')
dev.off()
####generazione a priori su k########
rzeta=function(a){
  b=2^{a-1}
  cond=TRUE
  while(cond){
    u=runif(1)
    v=runif(1)
    x=floor(u^(-1/(a-1)))
    if (x==Inf) x=.Machine$double.xmax
    t=(1+1/x)^(a-1)
    cond=(v*x*(t-1)/(b-1)>t/b)}
  return(x)}

NN2=c()
for(i in 1:100000) NN2[i]=rzeta(1.02)

k2=c()
for(i in 1:100000) {
  if (NN2[i]>.Machine$integer){
    k2[i]=length(unique(sample(.Machine$integer,size=500,rep=T)))
  }else{
    k2[i]=length(unique(sample(NN2[i],size=500,rep=T)))
  }
}

pdf("reFigure2.pdf")
par(mfrow = c(1,2))

#posteriori di k
NZ <- result$k1_chain[-(1:500)]
p.post=table(NZ)/length(NZ)
x.post=as.numeric(names(table(NZ)))
plot(x.post,p.post,xlim=c(410,495),xlab="k",type="b",lwd=1,cex.lab=2,main="",ylab="")
p.prior=table(k2)/length(k2)
x.prior=as.numeric(names(table(k2)))
lines(x.prior,p.prior,col=2,lwd=2)
hist(result$N1_chain)
dev.off()

getCouplingTime <- function(nRepeats = 1000){
  mTimes <- rep(NA, nRepeats)
  for (i in 1:nRepeats){
    mTimes[i] <- coupleGibbs(nMCMC, lambda1, lambda2, N1, N2, g, p, a1, a2,n, earlyStop = TRUE)$meetTime
  }
  return(mTimes)
}
mTimes <- getCouplingTime(nRepeats = nRepeats)
pdf("meetingTime.pdf")
hist(mTimes, "histogram of coupling time")
dev.off()
