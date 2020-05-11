### this function performs a full sweep update for theta 
### it updates each field of theta independently


update_theta <- function(theta, partition, partition_ll, alpha){
  theta_accept <- rep(0, length(theta))
  ## for each field
  for (l in 1:p){
    ## current state
    x <- theta[[l]] 
    cl_log_lik <- partition_ll[,l] 
    ## dirichlet proposal
    x_propose <- gtools::rdirichlet(1, alpha = (1 + concentration * x))[1,]
    cl_log_lik_new <- compute_loglikelihood_all_clusters_one_field_cpp(l - 1, partition, x_propose, V - 1, alpha)
    ## transition ratio 
    lratio <- log(gtools::ddirichlet(x = x, alpha = (1 + concentration * x_propose))) - log(gtools::ddirichlet(x = x_propose, alpha = (1 + concentration * x)))
    ## log accept probability
    laccept <- lratio - sum(cl_log_lik[which(partition$clsize != 0)]) + sum(cl_log_lik_new[which(partition$clsize != 0)])
    if (log(runif(1)) < laccept){
      theta[[l]] <- x_propose
      theta_accept[l] <- 1
      partition_ll[,l] <- cl_log_lik_new
    }
  }
  return(list(theta = theta, theta_accept = theta_accept, partition_ll = partition_ll))
}
