## sample from a maximum coupling of Normal(mu1, sigma) and Normal(mu2, sigma) 
## where mu1, mu2 are the means and sigma is the standard deviation
## The function returns a list where 
## 'xy' contains the two draws in a vector c(x,y)
## 'identical' is TRUE if x == y
rnorm_reflectionmxax <- function (mu1, mu2, sigma) 
{
  reflmax_xy <- c(0, 0)
  xdot <- rnorm(1)
  z <- (mu1 - mu2)/sigma
  normz <- sqrt(sum(z^2))
  e <- z/normz
  utilde <- runif(1, 0, 1)
  accept <- (log(utilde) < (dnorm(xdot + z, 0, 1, log = TRUE) - 
                              dnorm(xdot, log = TRUE)))
  ident_ <- FALSE
  if (accept) {
    ydot <- xdot + z
    ident_ <- TRUE
  }
  else {
    ydot <- xdot - 2 * (e * xdot) * e
  }
  reflmax_xy[1] <- mu1 + sigma * xdot
  reflmax_xy[2] <- mu2 + sigma * ydot
  return(list(xy = reflmax_xy, identical = ident_))
}

# the two functions below give sampling from the maximal coupling where the standard deviations 
# are not necessarily the same

get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)){
      return(list(xy = c(x,x), identical = TRUE))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq(1)
        reject <- (dq(y) + log(runif(1)) < dp(y))
      }
      return(list(xy = c(x,y), identical = FALSE))
    }
  }
}

rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
                        function(x) dnorm(x, mu1, sigma1, log = TRUE),
                        function(n) rnorm(n, mu2, sigma2),
                        function(x) dnorm(x, mu2, sigma2, log = TRUE))
  return(f())
}