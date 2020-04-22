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

# mu1 <- 1
# mu2 <- 2
# sigma <- 0.5
# 
# coupledpoints <- lapply(1:1e4, function(index) rnorm_reflectionmxax(mu1, mu2, sigma))
# mean(sapply(coupledpoints, function(l) l$identical))
# coupledpoints_x <- sapply(coupledpoints, function(l) l$xy[1])
# coupledpoints_y <- sapply(coupledpoints, function(l) l$xy[2])
# 
# hist(coupledpoints_x, prob = TRUE, nclass = 50, xlim = c(-2,5), main = '', xlab = '')
# curve(dnorm(x, mu1, sigma), add = TRUE)
# hist(coupledpoints_y, prob = TRUE, nclass = 50, add = TRUE, col = rgb(1,0,0,0.5))
# curve(dnorm(x, mu2, sigma), add = TRUE, col = rgb(1,0,0))
