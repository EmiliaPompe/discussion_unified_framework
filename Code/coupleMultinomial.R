## maximal coupling of two multinomial distributions
coupleMultinomial <- function(p, q, n){
	## p and q are a probability vectors of length n 
	## returns a sample from the maximal coupling of 
	## multinomial distributions specified by p and q
	x <- sample.int(n = n, size = 1, prob = p)
	w <- runif(1)
	if (w < q[x] / p[x]){
	  return(c(x,x))
	}else{
	  accept <- FALSE
	  while(! accept){
	    y <- sample.int(n = n, size = 1, prob = q)
	    w <- runif(1)
	    accept <- (w > p[y] / q[y])
	  }
	  return(c(x,y))
	}
}