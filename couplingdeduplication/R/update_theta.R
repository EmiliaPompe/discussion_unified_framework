### this function performs a full sweep update for theta 
### it updates each field of theta independently

# #'@export
#' update_theta <- function(theta, partition, partition_ll, alpha, concentration){
#'   theta_accept <- rep(0, length(theta))
#'   ## for each field
#'   for (l in 1:p){
#'     ## current state
#'     x <- theta[[l]] 
#'     cl_log_lik <- partition_ll[,l] 
#'     ## dirichlet proposal
#'     x_propose <- gtools::rdirichlet(1, alpha = (1 + concentration * x))[1,]
#'     logx_propose <- log(x_propose)
#'     cl_log_lik_new <- compute_loglikelihood_all_clusters_one_field_cpp(l - 1, partition, x_propose, logx_propose, V - 1, alpha)
#'     ## transition ratio 
#'     lratio <- log(gtools::ddirichlet(x = x, alpha = (1 + concentration * x_propose))) - log(gtools::ddirichlet(x = x_propose, alpha = (1 + concentration * x)))
#'     ## log accept probability
#'     laccept <- lratio - sum(cl_log_lik[which(partition$clsize != 0)]) + sum(cl_log_lik_new[which(partition$clsize != 0)])
#'     if (log(runif(1)) < laccept){
#'       theta[[l]] <- x_propose
#'       theta_accept[l] <- 1
#'       partition_ll[,l] <- cl_log_lik_new
#'     }
#'   }
#'   return(list(theta = theta, theta_accept = theta_accept, partition_ll = partition_ll))
#' }
# 
# 
# x <- c(A = 0.086, B = 0.042, C = 0.036, D = 0.028, E = 0.036, 
#            F = 0.032, G = 0.088, H = 0.108, I = 0.036, J = 0.038, K = 0.058, 
#            L = 0.002, M = 0.118, N = 0.008, O = 0.004, P = 0.034, R = 0.032, 
#            S = 0.104, T = 0.026, U = 0.038, V = 0.002, W = 0.044)
# concentration <- 1e4
# x_propose <- gtools::rdirichlet(1, alpha = (1 + concentration * x))[1,]
# 
# x_other <- rgamma(n = length(x), shape = concentration*x+1, rate = 1)
# x_other/sum(x_other)
# 
# 
# log(gtools::ddirichlet(x = x, alpha = (1 + concentration * x_propose))) - log(gtools::ddirichlet(x = x_propose, alpha = (1 + concentration * x)))
# lqratio <- 0
# sum((concentration*x) * log(x_propose) - (concentration*x_propose) * log(x)) + sum(lgamma(concentration*x_propose+1) - lgamma(concentration*x+1))

