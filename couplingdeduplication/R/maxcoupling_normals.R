#'@rdname maxcoupling_normal
#'@title Maximal coupling of two univariate Normal distributions
#'@description Sample from maximal coupling of two univariate Normal distributions,
#'specified through their means and standard deviations. See \code{\link{rmvnorm_max}} for a multivariate version.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param sigma1 First mean
#'@param sigma2 Second mean
#'@return Returns a list with
#'
#' \itemize{
#' \item "xy": the pair of samples \eqn{(x,y)}
#'
#' \item "identical": TRUE if \eqn{x = y}, FALSE otherwise
#' }
#'@export
maxcoupling_normal <- function(mu1, mu2, sigma1, sigma2){
  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
                        function(x) dnorm(x, mu1, sigma1, log = TRUE),
                        function(n) rnorm(n, mu2, sigma2),
                        function(x) dnorm(x, mu2, sigma2, log = TRUE))
  return(f())
}

# reflection-maximal coupling in one dimension
#'@rdname reflectionmaxcoupling
#'@title Reflection-maximal coupling of two univariate Normal distributions
#'@description Sample from reflection-maximal coupling of two univariate Normal distributions,
#'specified through their means, with common standard deviation.
#'See \code{\link{rmvnorm_reflectionmax}} for a multivariate version.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param sigma Common standard deviation
#'@return Returns a list with
#' \itemize{
#'
#' \item "xy": the pair of samples \eqn{(x,y)}
#'
#' \item "identical": TRUE if \eqn{x = y}, FALSE otherwise
#' }
#'@export
reflectionmaxcoupling <- function(mu1, mu2, sigma){
  # number of samples
  reflmax_xy <- c(0,0)
  # draw std normal first
  xdot <- rnorm(1)
  # this follows the notation of Bou Rabee et al, 2018, roughly
  z <- (mu1 - mu2) / sigma
  normz <- sqrt(sum(z^2))
  e <- z / normz
  utilde <- runif(1, 0, 1)
  accept <- (log(utilde) < (dnorm(xdot + z, 0, 1, log = TRUE) - dnorm(xdot, log = TRUE)))
  ident_ <- FALSE
  if (accept){
    ydot <- xdot + z
    ident_ <- TRUE
  } else {
    ydot <- xdot - 2 * (e * xdot) * e
  }
  reflmax_xy[1] <- mu1 + sigma * xdot
  reflmax_xy[2] <- mu2 + sigma * ydot
  return(list(xy = reflmax_xy, identical = ident_))
}


