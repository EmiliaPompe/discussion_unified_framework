

couple_multinomial_alt <- function(p, q, n){
  ## p and q are a probability vectors of length n 
  ## returns a sample from the maximal coupling of 
  ## multinomial distributions specified by p and q
  minpq <- pmin(p,q)
  sum_minpq <- sum(minpq)
  #
  if (runif(1) < sum_minpq){
    x <- sample.int(n = n, size = 1, prob = minpq / sum_minpq)
    return(c(x,x))
  } else {
    x <- sample.int(n = n, size = 1, prob = (p - minpq)/(1 - sum_minpq))
    y <- sample.int(n = n, size = 1, prob = (q - minpq)/(1 - sum_minpq))
    return(c(x,y))
  }
}

# n <- 14
#
# p <- rexp(n, 1) 
# p <- p/sum(p)
# q <- rexp(n, 1)
# q <- q/sum(q)
#
# xy <- t(sapply(1:1e5, function(rep) couple_multinomial_alt(p,q,n)))
# summary(as.numeric((table(xy[,1 ])/nrow(xy) - p)/p))
# summary(as.numeric((table(xy[,2 ])/nrow(xy) - q)/q))
