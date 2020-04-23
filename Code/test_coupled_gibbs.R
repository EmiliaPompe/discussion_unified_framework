#### This file reproduces the two panels of Figure 2 in the paper
#### coupled updates of lambda and N and theta
#### keeping alpha fixed
rm(list = ls())
set.seed(2020)
source("getData.R")
source("initPara.R")
library(Rcpp)
sourceCpp("init_clustering.cpp")
sourceCpp("compute_loglikelihood_clusters.cpp")
sourceCpp("compute_psampq.cpp")
sourceCpp("computeNweights.cpp")
source("couple_dirichlet.R")
source("couple_multinomial_alt.R")
source("couple_normal.R")
source("coupleTheta.R")
source("coupleLambda.R")
source("coupleN.R")
source("relabel.R")
source("coupledGibbs.R")

result <- coupleGibbs(nMCMC, N1, N2, g, p1,p2, a1, a2,n, earlyStop = FALSE)

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
    result <- coupleGibbs(nMCMC, N1, N2, g, p1, p2, a1, a2,n, earlyStop = TRUE)
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

# pdf("distance.pdf", width = 6, height = 4)
# par(mfrow = c(1,3), mar = c(4,4,1,1))
# title <- paste( "histogram of coupling time, with", nRepeats, "repeats", sep = " ")
# hist(repeatCoupling$meeting_times, main = title, xlab = 'meetimg times')
# plot(ts, tv_upper, main = '', type = 'l',
#      xlab = paste('iteration, lag =', L, sep = ' '),
#      ylab = 'upper bound d_TV')
# abline(h = 1,col = 'red')
# plot(ts, colMeans(w_upper[,ts]),
#      xlab = paste('iteration, lag =', L, sep = ' '),
#      ylab = 'upper bound d_W',
#      type = 'l')
# dev.off()

