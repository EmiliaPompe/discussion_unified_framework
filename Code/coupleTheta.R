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
    x_propose <- gtools::rdirichlet(1, alpha = concentraton * x)[1,]
    p_propose[ cumdime[l] : (cumdime[l+1] - 1)] <- x_propose
    ## transition ratio 
    lratio <- log(gtools::ddirichlet(x = x, alpha = concentraton * x_propose)) - log(gtools::ddirichlet(x = x_propose, alpha = concentraton * x))
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

coupleTheta <- function(p1, clustering1,
                        p2, clustering2){
  cl_log_lik1 <- compute_loglikelihood_clusters(clustering1, p1, V, dimV, a)
  p1_propose <- p1
  cl_log_lik2 <- compute_loglikelihood_clusters(clustering2, p2, V, dimV, a)
  p2_propose <- p2
  for ( l in 1:H){
    ## first propose from coupled dirichlet
    x1 <- p1[cumdime[l] : (cumdime[l+1] - 1)] 
    x2 <- p1[cumdime[l] : (cumdime[l+1] - 1)]
    cp_dirichlet <- couple_dirichlet(mean1 = x1, mean2 = x2, inversescale = 1 / concentraton)$xy
    x1_propose <- cp_dirichlet[1:length(x1)]
    x2_propose <- cp_dirichlet[(length(x1) + 1) : (2 * length(x1))]
    p1_propose[ cumdime[l] : (cumdime[l+1] - 1)] <- x1_propose
    p2_propose[ cumdime[l] : (cumdime[l+1] - 1)] <- x2_propose
    lratio1 <- log(gtools::ddirichlet(x = x1, alpha = concentraton * x1_propose)) - log(gtools::ddirichlet(x = x1_propose, alpha = concentraton * x1))
    lratio2 <- log(gtools::ddirichlet(x = x2, alpha = concentraton * x2_propose)) - log(gtools::ddirichlet(x = x2_propose, alpha = concentraton * x2))
    cl_log_lik1_new <- compute_loglikelihood_clusters(clustering1, p1_propose, V, dimV, a)
    cl_log_lik2_new <- compute_loglikelihood_clusters(clustering2, p2_propose, V, dimV, a)
    laccept1 <- lratio1 - sum_notInf(cl_log_lik1[,l]) + sum_notInf(cl_log_lik1_new[,l])
    laccept2 <- lratio2 - sum_notInf(cl_log_lik2[,l]) + sum_notInf(cl_log_lik2_new[,l])
    lU <- log(runif(1))
    if (lU < laccept1){
      p1 <- p1_propose
      cl_log_lik1 <- cl_log_lik1_new
    }
    if (lU < laccept2){
      p2 <- p2_propose
      cl_log_lik2 <- cl_log_lik2_new
    }
  }
  return(list(p1 = p1, p2 = p2))
}

# source("getData.R")
# source("initPara.R")
# source("couple_dirichlet.R")
# library(Rcpp)
# sourceCpp("init_clustering.cpp")
# ## lambda is eta  in the paper
# lambda <- sample(n, size = n, replace = T)
# print(lambda)
# clustering <- init_clustering(lambda-1)
# 
# cl_log_lik <- compute_loglikelihood_clusters(clustering, p, V, dimV, a)
# 
# # update_theta(p, clustering)
# 
# # sum(update_theta(p, clustering)[1:22])
# lambda1 <- sample(n, size = n, replace = T)
# lambda2 <- sample(n, size = n, replace = T)
# p1 <- p
# p2 <- rep(NA, length(p1))
# for (l in 1:H){
#   p2[cumdime[l] : (cumdime[l+1] -1 ) ] <- gtools::rdirichlet(1, alpha = p1[cumdime[l] : (cumdime[l+1] - 1)] * concentraton )[1,]
#   # print(sum(p2[cumdime[l] : (cumdime[l+1] -1 ) ]))
# }
# clustering1 <- init_clustering(lambda1 -1 )
# clustering2 <- init_clustering(lambda2 - 1)
# 
 
# result <- coupleTheta(p1, clustering1,
#                         p2, clustering2)
# result
