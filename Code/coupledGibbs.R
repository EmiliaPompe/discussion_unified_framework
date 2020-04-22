#### This file reproduces the two panels of Figure 2 in the paper
#### coupled updates of lambda and N, keeping  theta fixed
#### it does not update alpha, but alpha changes in each iteration due to re-labeling
#### author: Phyllis Ju
#### email: nju@g.harvard.edu
#### updated: 21-April 2020

rm(list = ls())
source("getData.R")
source("initPara.R")
# source("coupleMultinomial.R")
library(Rcpp)
sourceCpp("computeNweights.cpp")
# sourceCpp("computeQcpp.cpp")
sourceCpp("init_clustering.cpp")
sourceCpp("compute_psampq.cpp")
source("coupleLambda.R")
source("couple_multinomial_alt.R")
source("coupleN.R")
source("relabel.R")


coupleGibbs <- function(nMCMC, N1, N2, g, p, a1, a2, n, earlyStop = TRUE){
  lambda1 <- sample(n, size = n, replace = T)
  lambda2 <- sample(n, size = n, replace = T)
  lambda1 <- relabel(lambda1)$lambda
  lambda2 <- relabel(lambda2)$lambda
  meetTime <- Inf
  N1_chain <- N2_chain <- rep(0, nMCMC)
  k1_chain <- k2_chain <- rep(0, nMCMC)
  lambda1_chain <- lambda2_chain <- matrix(0, nrow = n, ncol = nMCMC)
  percentageCoupled <- rep(0, nMCMC)
  w_upper <- rep(0,nMCMC)
  L <- 1
  ### run L iterations for the chain lambda1, N1 and a1
  for (iter in 1:L){
    for (j in 1:n){
      lambdajs <- coupleLambda(j, lambda1, lambda2, p1 = p, p2 = p, a1 = a1, a2 = a2, N1 = N1, N2 = N2)
      lambda1[j] <- lambdajs[1]
    }
    relabel1 <- relabel(lambda1)
    lambda1 <- relabel1$lambda
    a1 <- a1[relabel1$iis,]
    ksize1 <- length(unique(lambda1))
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize1, g, n)
    N1 <- Ns[1]
    ### store 
    k1_chain[iter] <- ksize1
    N1_chain[iter] <- N1
    lambda1_chain[, iter] <- lambda1
  }
  ### try to couple the chains with lag L
  iter <- L
  for (iter in (L + 1):nMCMC){
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
    # ## relabel lambd1 and lambda2 and change a1 a2 accordingly
    relabel1 <- relabel(lambda1)
    lambda1 <- relabel1$lambda
    a1 <- a1[relabel1$iis,]
    relabel2 <- relabel(lambda2)
    lambda2 <- relabel2$lambda
    a2 <- a2[relabel2$iis,]
    lambda1_chain[, iter] <- lambda1
    lambda2_chain[, iter - L] <- lambda2 
    if (iter %% (nMCMC / 50) == 0 ){
      print(paste('iteration', iter, percentageCoupled[iter] * 100, "percent coupled",sep = ' '))
    }
    # partition1 <- init_clustering(lambda1 - 1)
    # partition2 <- init_clustering(lambda2 - 1)
    ksize1 <- length(unique(lambda1))
    ksize2 <- length(unique(lambda2))
    k1_chain[iter] <- ksize1
    k2_chain[iter - L] <- ksize2
    ## update alpha 
    
    ## update theta
    # for (l in 1:H){
    #   
    # }
    ## update N 
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize2, g, n)
    N1 <- Ns[1]
    N2 <- Ns[2]
    N1_chain[iter] <- N1
    N2_chain[iter - L] <- N2
    meet <- FALSE
    if (N1 == N2 & all(lambda1 == lambda2) & !meet){
      meetTime <- iter
      meet <- TRUE
    }
    if (!meet & earlyStop){
      w_upper[iter - L] <- abs(N1_chain[iter] - N2_chain[iter - L]) + sum(abs(lambda1_chain[, iter] - lambda2_chain[, iter - L]))
    }
    if (earlyStop & meet){
      break
    }
  }
  return(list(N1_chain = N1_chain, N2_chain = N2_chain,
              k1_chain = k1_chain, k2_chain = k2_chain,
              meetTime = meetTime, percentageCoupled = percentageCoupled,
              lambda1_chain = lambda1_chain, lambda2_chain = lambda2_chain,
              L = L, w_upper = w_upper))
}

result <- coupleGibbs(nMCMC, N1, N2, g, p, a1, a2,n, earlyStop = FALSE)
# print(result$meetTime)
# print(which(result$N1_chain != result$N2_chain))
# 
# 
# # pdf("reFigure2.pdf")
# par(mfrow = c(1,2))
# ## save this plot
# hist(result$k1_chain)
# hist(result$N1_chain)
# # dev.off()
# par(mfrow = c(1,2))
# plot(result$N1_chain, type = 'l')
# lines(result$N2_chain, col = 'red')
# plot(result$k1_chain, type = 'l')
# lines(result$k2_chain, col = 'red')
# 
# pdf("percentageCoupledLambda.pdf")
# par(mfrow = c(1,1))
# plot(result$percentageCoupled, type  = 'l')
# dev.off()
# ####generazione a priori su k########
# rzeta=function(a){
#   b=2^{a-1}
#   cond=TRUE
#   while(cond){
#     u=runif(1)
#     v=runif(1)
#     x=floor(u^(-1/(a-1)))
#     if (x==Inf) x=.Machine$double.xmax
#     t=(1+1/x)^(a-1)
#     cond=(v*x*(t-1)/(b-1)>t/b)}
#   return(x)}
# 
# NN2=c()
# for(i in 1:100000) NN2[i]=rzeta(1.02)
# 
# k2=c()
# for(i in 1:100000) {
#   if (NN2[i]>.Machine$integer){
#     k2[i]=length(unique(sample(.Machine$integer,size=500,rep=T)))
#   }else{
#     k2[i]=length(unique(sample(NN2[i],size=500,rep=T)))
#   }
# }

# pdf("reFigure2.pdf")
# par(mfrow = c(1,2))
# 
# #posteriori di k
# NZ <- result$k1_chain[-(1:nBurn)]
# p.post=table(NZ)/length(NZ)
# x.post=as.numeric(names(table(NZ)))
# plot(x.post,p.post,xlim=c(410,495),xlab="k",type="b",lwd=1,cex.lab=2,main="",ylab="")
# p.prior=table(k2)/length(k2)
# x.prior=as.numeric(names(table(k2)))
# lines(x.prior,p.prior,col=2,lwd=2)
# hist(result$N1_chain)
# dev.off()

getCouplingTime <- function(nRepeats = 1000){
  mTimes <- rep(NA, nRepeats)
  wUpper <- matrix(data = 0, nrow = nRepeats, ncol = nMCMC)
  for (i in 1:nRepeats){
    result <- coupleGibbs(nMCMC, N1, N2, g, p, a1, a2,n, earlyStop = TRUE)
    mTimes[i] <- min(nMCMC, result$meetTime)
    ## this is equation 4 in Biswas et al. 2019
    for (t in 1 : (mTimes[i] - 1 - result$L)){
      wUpper[i, t] <- sum(result$w_upper[1: (mTimes[i] - result$L - t)])
    }
  }
  return(list(meeting_times = mTimes,
              wUpper = wUpper))
}

repeatCoupling <- getCouplingTime(nRepeats = nRepeats)
# pdf("meetingTime.pdf")
title <- paste( "histogram of coupling time, with", nRepeats, "repeats", sep = " ")
hist(repeatCoupling$meeting_times, main = title, xlab = 'meetimg times')
# dev.off()

print(table(repeatCoupling$meeting_times))

## process the tv_upper bounds 
max_meeting <- max(repeatCoupling$meeting_times)
ts <- c(1: min(nMCMC, (2 + max_meeting)))
tv_upper <- rep(NA, length(ts))
for (t in ts){
  tv_upper[t ] <- mean( pmax(0, repeatCoupling$meeting_times - t))
}
par(mfrow = c(1,2), mar = c(4,4,1,1))
plot(ts, tv_upper, main = '', type = 'l',
     xlab = 'iteration',
     ylab = 'upper bound d_TV',
     log = 'x')
abline(h = 1,col = 'red')
w_upper <- repeatCoupling$wUpper
plot(colMeans(w_upper[,ts]),
     xlab = paste('iteration, lag =', L, sep = ' '),
     ylab = 'upper bound d_W',
     type = 'l')

pdf("distance.pdf", width = 6, height = 4)
par(mfrow = c(1,3), mar = c(4,4,1,1))
title <- paste( "histogram of coupling time, with", nRepeats, "repeats", sep = " ")
hist(repeatCoupling$meeting_times, main = title, xlab = 'meetimg times')
plot(ts, tv_upper, main = '', type = 'l',
     xlab = paste('iteration, lag =', L, sep = ' '),
     ylab = 'upper bound d_TV')
abline(h = 1,col = 'red')
plot(ts, colMeans(w_upper[,ts]),
     xlab = paste('iteration, lag =', L, sep = ' '),
     ylab = 'upper bound d_W',
     type = 'l')
dev.off()
