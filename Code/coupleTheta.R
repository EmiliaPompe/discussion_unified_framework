sum_notInf <- function(xx){
	sum(xx[which(is.finite(xx))])
}

update_theta <- function(p, clustering){
  ## for each field
	cl_log_lik <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)
	p_propose <- p
  for (l in 1:H){
  	## current state
    x <- p[cumdime[l] : (cumdime[l+1] - 1)] 
    ## dirichlet proposal
    x_propose <- gtools::rdirichlet(1, alpha = (1 + concentration * x))[1,]
    p_propose[ cumdime[l] : (cumdime[l+1] - 1)] <- x_propose
    ## transition ratio 
    lratio <- log(gtools::ddirichlet(x = x, alpha = (1 + concentratoin * x_propose))) - log(gtools::ddirichlet(x = x_propose, alpha = (1 + concentration * x)))
    cl_log_lik_new <- compute_loglikelihood_clusters(clustering, p_propose, V, dimV, a)
    ## log accept probability
    laccept <- lratio - sum_notInf(cl_log_lik[,l]) + sum_notInf(cl_log_lik_new[,l])
    if (log(runif(1)) < laccept){
    	p <- p_propose
    	cl_log_lik <- cl_log_lik_new
    }
  }
  return(p)
}

# clustering1 <- partition1
# clustering2 <- partition2

coupleTheta <- function(theta1, theta2, partition1, partition2, alpha1, alpha2){
  for ( l in 1 : L){
    ## first propose from coupled dirichlet
    x1 <- theta1[cumdime[l] : (cumdime[l+1] - 1)] 
    x2 <- theta2[cumdime[l] : (cumdime[l+1] - 1)]
    cp_dirichlet <- couple_dirichlet(mean1 = (1 / concentration + x1) , mean2 = (1 / concentration + x2), inversescale = 1 / concentration)$xy
    x1_propose <- cp_dirichlet[1:length(x1)]
    x2_propose <- cp_dirichlet[(length(x1) + 1) : (2 * length(x1))]
    # print(all(x1_propose == x2_propose))
    lratio1 <- log(gtools::ddirichlet(x = x1, alpha = (1 + concentration * x1_propose))) - log(gtools::ddirichlet(x = x1_propose, alpha = (1 + concentration * x1)))
    lratio2 <- log(gtools::ddirichlet(x = x2, alpha = (1 + concentration * x2_propose))) - log(gtools::ddirichlet(x = x2_propose, alpha = (1 + concentration * x2)))
    cl_field_loglik1 <- compute_loglikelihood_clusters_field(l - 1, partition1, x1, V, dimV, alpha1)
    cl_field_loglik2 <- compute_loglikelihood_clusters_field(l - 1, partition2, x2, V, dimV, alpha2)
    cl_field_loglik1_propose <- compute_loglikelihood_clusters_field(l - 1, partition1, x1_propose, V, dimV, alpha1)
    cl_field_loglik2_propose <- compute_loglikelihood_clusters_field(l - 1, partition2, x2_propose, V, dimV, alpha2)
    laccept1 <- lratio1 - sum_notInf(cl_field_loglik1) + sum_notInf(cl_field_loglik1_propose)
    laccept2 <- lratio2 - sum_notInf(cl_field_loglik2) + sum_notInf(cl_field_loglik2_propose)
    lU <- log(runif(1))
    if (lU < laccept1){
      theta1[ cumdime[l] : (cumdime[l + 1] -1)]  <- x1_propose
    }
    if (lU < laccept2){
      theta2[ cumdime[l] : (cumdime[l + 1] -1)]  <- x2_propose
    }
  }
  return(list(theta1 = theta1, theta2 = theta2))
}

# test <- T
# source("getData.R")
# source("initPara.R")
# # source("coupleMultinomial.R")
# library(Rcpp)
# # sourceCpp("computeNweights.cpp")
# # sourceCpp("computeQcpp.cpp")
# # source("coupleLambda.R")
# # source("couple_multinomial_alt.R")
# # source("coupleN.R")
# # source("relabel.R")
# sourceCpp("init_clustering.cpp")
# sourceCpp("compute_loglikelihood_clusters.cpp")
# sourceCpp('compute_loglikelihood_clusters_field.cpp')
# source("couple_dirichlet.R")
# 
# for ( l in 1 : L){
#   ## first propose from coupled dirichlet
#   x1 <- theta1[cumdime[l] : (cumdime[l+1] - 1)] 
#   x2 <- theta2[cumdime[l] : (cumdime[l+1] - 1)]
#   cp_dirichlet <- couple_dirichlet(mean1 = (1 / concentration + x1) , mean2 = (1 / concentration + x2), inversescale = 1 / concentration)$xy
#   x1_propose <- cp_dirichlet[1:length(x1)]
#   x2_propose <- cp_dirichlet[(length(x1) + 1) : (2 * length(x1))]
#   print(all(x1_propose == x2_propose))
#   print(x1_propose)
#   print(x2_propose)
#   lratio1 <- log(gtools::ddirichlet(x = x1, alpha = (1 + concentration * x1_propose))) - log(gtools::ddirichlet(x = x1_propose, alpha = (1 + concentration * x1)))
#   lratio2 <- log(gtools::ddirichlet(x = x2, alpha = (1 + concentration * x2_propose))) - log(gtools::ddirichlet(x = x2_propose, alpha = (1 + concentration * x2)))
#   cl_field_loglik1 <- compute_loglikelihood_clusters_field(l - 1, partition1, x1, V, dimV, alpha1)
#   cl_field_loglik2 <- compute_loglikelihood_clusters_field(l - 1, partition2, x2, V, dimV, alpha2)
#   cl_field_loglik1_propose <- compute_loglikelihood_clusters_field(l - 1, partition1, x1_propose, V, dimV, alpha1)
#   cl_field_loglik2_propose <- compute_loglikelihood_clusters_field(l - 1, partition2, x2_propose, V, dimV, alpha2)
#   # print(cl_field_loglik1)
#   # print(compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)[,l])
#   # print(cl_field_loglik2)
#   # print(compute_loglikelihood_clusters(partition2, theta2, V, dimV, alpha2)[,l])
#   laccept1 <- lratio1 - sum_notInf(cl_field_loglik1) + sum_notInf(cl_field_loglik1_propose)
#   laccept2 <- lratio2 - sum_notInf(cl_field_loglik2) + sum_notInf(cl_field_loglik2_propose)
#   lU <- log(runif(1))
#   if (lU < laccept1){
#     theta1[ cumdime[l] : (cumdime[l + 1] -1)]  <- x1_propose
#   }
#   if (lU < laccept2){
#     theta2[ cumdime[l] : (cumdime[l + 1] -1)]  <- x2_propose
#   }
# }
# 
