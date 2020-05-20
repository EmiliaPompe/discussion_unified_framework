coupleGibbs <- function(nMCMC, N1, N2, g, theta1, theta2, alpha1, alpha2, n, lag, earlyStop = TRUE){
  eta1 <- sample(n, size = n, replace = T)
  eta2 <- sample(n, size = n, replace = T)
  eta1 <- relabel(eta1)$lambda
  eta2 <- relabel(eta2)$lambda
  partition1 <- init_clustering(eta1 - 1)
  partition2 <- init_clustering(eta2 - 1)
  cl_log_lik1 <- compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)
  cl_log_lik2 <- compute_loglikelihood_clusters(partition2, theta2, V, dimV, alpha2)
  meetTime <- Inf
  N1_chain <- N2_chain <- rep(0, nMCMC)
  k1_chain <- k2_chain <- rep(0, nMCMC)
  theta1_chain <- theta2_chain <- matrix(0, nrow = length(theta1), ncol = nMCMC)
  eta1_chain <- eta2_chain <- matrix(0, nrow = n, ncol = nMCMC)
  percentageCoupled <- rep(0, nMCMC)
  w_upper <- rep(0,nMCMC)
  ### run lag iterations for the chain lambda1, N1 and a1
  for (iter in 1:lag){
    updateLambdaRes <- update_lambda(eta1 - 1, partition1, cl_log_lik1, theta1, V, dimV, alpha1, N1)
    eta1 <- updateLambdaRes$lambda + 1
    partition1 <- updateLambdaRes$clustering
    cl_log_lik1  <- updateLambdaRes$clusterloglikelihoods
    ksize1 <- length(unique(eta1))
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize1, g, n)
    N1 <- Ns[1]
    ### store 
    k1_chain[iter] <- ksize1
    N1_chain[iter] <- N1
    eta1_chain[, iter] <- eta1
  }
  cl_log_lik1 <- compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)
  ### try to couple the chains with lag L
  iter <- lag + 1
  meet <- FALSE
  for (iter in (lag + 1):nMCMC){
    ## update lambda 
    for (j in 1:n){
      cp_lambda_result <- couple_eta(j, eta1, eta2, partition1, partition2, cl_log_lik1, cl_log_lik2, theta1, theta2, alpha1, alpha2, N1, N2)
      eta1 <- cp_lambda_result$eta1
      eta2 <- cp_lambda_result$eta2
      partition1 <- cp_lambda_result$partition1
      partition2 <- cp_lambda_result$partition2
      cl_log_lik1 <- cp_lambda_result$cl_log_lik1
      cl_log_lik2 <- cp_lambda_result$cl_log_lik2
    }
    ## relabel lambd1 and lambda2 and change a1 a2 accordingly
    relabel1 <- relabel(eta1)
    eta1 <- relabel1$lambda
    alpha1 <- alpha1[relabel1$iis,]
    relabel2 <- relabel(eta2)
    eta2 <- relabel2$lambda
    alpha2 <- alpha2[relabel2$iis,]
    eta1_chain[, iter] <- eta1
    eta2_chain[, iter - lag] <- eta2 
    partition1 <- init_clustering(eta1 - 1)
    partition2 <- init_clustering(eta2 - 1)
    cl_log_lik1 <- compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)
    cl_log_lik2 <- compute_loglikelihood_clusters(partition2, theta2, V, dimV, alpha2)
    ksize1 <- length(unique(eta1))
    ksize2 <- length(unique(eta2))
    k1_chain[iter] <- ksize1
    k2_chain[iter - lag] <- ksize2
    ## update alpha 
    ## update theta
    # couple_theta <- coupleTheta(theta1, theta2,  partition1, partition2, alpha1, alpha2)
    # theta1 <- couple_theta$theta1
    # theta2 <- couple_theta$theta2
    # if(any(theta1 == 0) | any(theta2 == 0)){
    #   print(paste('couple_theta gives 0 at iteration',iter,'\n',sep = " "))
    # }
    # cl_log_lik1 <- compute_loglikelihood_clusters(partition1, theta1, V, dimV, alpha1)
    # cl_log_lik2 <- compute_loglikelihood_clusters(partition2, theta2, V, dimV, alpha2)
    theta1_chain[,iter] <- theta1
    theta2_chain[,iter - lag] <- theta2
    ## update N 
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize2, g, n)
    N1 <- Ns[1]
    N2 <- Ns[2]
    N1_chain[iter] <- N1
    N2_chain[iter - lag] <- N2
    percentageCoupled[iter] <- ( sum(eta1 == eta2) + sum(theta1 == theta2)) / (n + length(theta1))
    if (iter %% (nMCMC / 50) == 0 ){
      print(paste('iteration', iter, percentageCoupled[iter] * 100, "percent coupled",sep = ' '))
    }
    if (N1 == N2 & all(eta1 == eta2) & all(theta1 == theta2) & !meet){
      meetTime <- iter
      meet <- TRUE
    }
    if (!meet & earlyStop){
      w_upper[iter - lag] <- abs(N1_chain[iter] - N2_chain[iter - lag]) +
        sum(abs(eta1_chain[, iter] - eta2_chain[, iter - lag])) + 
        sum(abs(theta1 - theta2))
    }
    if (earlyStop & meet){
      break
    }
  }
  return(list(N1_chain = N1_chain, N2_chain = N2_chain,
              k1_chain = k1_chain, k2_chain = k2_chain,
              theta1_chain = theta1_chain, theta2_chain = theta2_chain, 
              meetTime = meetTime, percentageCoupled = percentageCoupled,
              eta1_chain = eta1_chain, eta2_chain = eta2_chain,
              lag = lag, w_upper = w_upper))
}

