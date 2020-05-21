## Script to try the multinomial sampler and the coupling of it 
## this is work in progress
rm(list = ls())
set.seed(2)
library(couplingdeduplication)

logw <- rnorm(8)
exp(logw)/sum(exp(logw))
nsamples <- 1e5
uniforms <- runif(nsamples)

draws <- sapply(1:nsamples, function(i) couplingdeduplication:::multinomial_(logw, uniforms[i]))
round(as.numeric(table(draws) / nsamples), 3)
round(exp(logw)/sum(exp(logw)), 3)


logw1 <- rnorm(8)
logw2 <- rnorm(8)
exp(logw)/sum(exp(logw))
uniforms1 <- runif(nsamples)
uniforms2 <- runif(nsamples)

draws <- t(sapply(1:nsamples, function(i) couplingdeduplication:::coupled_multinomial_(logw1, logw2, uniforms1[i], uniforms2[i])))
#
round(as.numeric(table(draws[,1]) / nsamples), 3)
round(exp(logw1)/sum(exp(logw1)), 3)
#
round(as.numeric(table(draws[,2]) / nsamples), 3)
round(exp(logw2)/sum(exp(logw2)), 3)

##
mean(apply(draws, 1, function(v) v[1]==v[2]))
sum(pmin(exp(logw1)/sum(exp(logw1)), exp(logw2)/sum(exp(logw2))))
