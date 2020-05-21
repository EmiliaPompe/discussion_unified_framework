## draw from max coupling of two Dirichlet with parameters alpha1 and alpha2
#'@export
maxcoupling_dirichlet <- function(alpha1, alpha2){
  x <- gtools::rdirichlet(1, alpha = alpha1)[1,]
  if (log(gtools::ddirichlet(x, alpha = alpha1)) + log(runif(1)) < log(gtools::ddirichlet(x, alpha = alpha2))){
    return(list(xy = cbind(x, x), identical = TRUE))
  }
  else {
    reject <- TRUE
    y <- NA
    while (reject) {
      y <- gtools::rdirichlet(1, alpha2)[1,]
      reject <- (log(gtools::ddirichlet(y, alpha = alpha2)) + log(runif(1)) < log(gtools::ddirichlet(y, alpha = alpha1)))
    }
    return(list(xy = cbind(x, y), identical = FALSE))
  }
}
