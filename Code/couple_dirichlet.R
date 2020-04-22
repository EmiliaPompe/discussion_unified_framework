## sample from a maximum coupling of two Dirichlet distributions
## centered at mean1 and mean2 (vectors in the simplex)
## and with inversescale parameter, which means these are Dirichlet 
## with 'alpha' = mean1 / inversescale, respectively mean2 / inversescale
## The smaller the 'inversecale', the more concentrated the draws are around their means.

## returns a list with a matrix 'xy', first column containing x and second y
## and 'identical' a boolean equals to TRUE if x == y

couple_dirichlet <- function (mean1, mean2, inversescale){
  x <- gtools::rdirichlet(1, alpha = mean1 / inversescale)[1,]
  if (log(gtools::ddirichlet(x, alpha = mean1 / inversescale)) + log(runif(1)) < log(gtools::ddirichlet(x, alpha = mean2 / inversescale))){
    return(list(xy = cbind(x, x), identical = TRUE))
  }
  else {
    reject <- TRUE
    y <- NA
    while (reject) {
      y <- gtools::rdirichlet(1, mean2 / inversescale)[1,]
      reject <- (log(gtools::ddirichlet(y, alpha = mean2 / inversescale)) + log(runif(1)) < log(gtools::ddirichlet(y, alpha = mean1 / inversescale)))
    }
    return(list(xy = cbind(x, y), identical = FALSE))
  }
}

# mean1 <- c(1/3,1/3,1/3)
# mean2 <- c(1/2,1/4,1/4)
# inversescale <- .04
# x1 <- gtools::rdirichlet(1e3, mean1 / inversescale)
# x2 <- gtools::rdirichlet(1e3, mean2 / inversescale)
# plot(x1[,1], x1[,2], xlim = c(0,1), ylim = c(0,1))
# points(x2[,1], x2[,2], col = 'red')
# 
# coupledpoints <- lapply(1:1e3, function(index) couple_dirichlet(mean1, mean2, inversescale))
# mean(sapply(coupledpoints, function(l) l$identical))
# coupledpoints_x <- t(sapply(coupledpoints, function(l) l$xy[,1]))
# coupledpoints_y <- t(sapply(coupledpoints, function(l) l$xy[,2]))
# plot(coupledpoints_x[,1], coupledpoints_x[,2], xlim = c(0,1), ylim = c(0,1))
# points(coupledpoints_y[,1], coupledpoints_y[,2], col = 'red')
