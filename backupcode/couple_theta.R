# sum_notInf <- function(xx){
# 	sum(xx[which(is.finite(xx))])
# }

couple_theta <- function(theta1, theta2, partition1, partition2, partition_ll1, partition_ll2,alpha1, alpha2){
  coupled <- TRUE
  for ( l in 1 : p){
    ## first propose from coupled dirichlet
    x1 <- theta1[[l]]
    x2 <- theta2[[l]]
    cp_dirichlet <- couple_dirichlet(mean1 = (1 / concentration + x1) , mean2 = (1 / concentration + x2), inversescale = 1 / concentration)$xy
    x1_propose <- cp_dirichlet[1:length(x1)]
    x2_propose <- cp_dirichlet[(length(x1) + 1) : (2 * length(x1))]
    lratio1 <- log(gtools::ddirichlet(x = x1, alpha = (1 + concentration * x1_propose))) - log(gtools::ddirichlet(x = x1_propose, alpha = (1 + concentration * x1)))
    lratio2 <- log(gtools::ddirichlet(x = x2, alpha = (1 + concentration * x2_propose))) - log(gtools::ddirichlet(x = x2_propose, alpha = (1 + concentration * x2)))
    cl_field_loglik1 <- partition_ll1[,l]
    cl_field_loglik2 <- partition_ll2[,l]
    cl_field_loglik1_propose <- compute_loglikelihood_all_clusters_one_field_cpp(l - 1, partition1, x1_propose, V - 1, alpha1)
    cl_field_loglik2_propose <- compute_loglikelihood_all_clusters_one_field_cpp(l - 1, partition2, x2_propose, V - 1, alpha2)
    laccept1 <- lratio1 - sum(cl_field_loglik1[which(partition1$clsize != 0)]) +sum(cl_field_loglik1_propose[which(partition1$clsize != 0)])
    laccept2 <- lratio2 - sum(cl_field_loglik2[which(partition2$clsize != 0)]) +sum(cl_field_loglik2_propose[which(partition2$clsize != 0)])
    lU <- log(runif(1))
    if (lU < laccept1){
      theta1[[l]] <- x1_propose
      partition_ll1[,l] <- cl_field_loglik1_propose
    }
    if (lU < laccept2){
      theta2[[l]]  <- x2_propose
      partition_ll2[,l] <- cl_field_loglik2_propose
    }
    coupled <- ( coupled & all(theta1[[l]] == theta2[[l]]))
  }
  return(list(theta1 = theta1, theta2 = theta2,
              partition_ll1 = partition_ll1, partition_ll2 = partition_ll2,
              coupled = coupled))
}
