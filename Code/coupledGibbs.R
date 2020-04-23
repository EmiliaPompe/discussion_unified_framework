coupleGibbs <- function(nMCMC, N1, N2, g, p1, p2, a1, a2, n, earlyStop = TRUE){
  lambda1 <- sample(n, size = n, replace = T)
  lambda2 <- sample(n, size = n, replace = T)
  lambda1 <- relabel(lambda1)$lambda
  lambda2 <- relabel(lambda2)$lambda
  partition1 <- init_clustering(lambda1 - 1)
  partition2 <- init_clustering(lambda2 - 1)
  ll_clusters1 <- compute_loglikelihood_clusters(partition1, p1, V, dimV, a1)
  ll_clusters2 <- compute_loglikelihood_clusters(partition2, p2, V, dimV, a2)
  meetTime <- Inf
  N1_chain <- N2_chain <- rep(0, nMCMC)
  k1_chain <- k2_chain <- rep(0, nMCMC)
  p1_chain <- p2_chain <- matrix(0, nrow = length(p1), ncol = nMCMC)
  lambda1_chain <- lambda2_chain <- matrix(0, nrow = n, ncol = nMCMC)
  percentageCoupled <- rep(0, nMCMC)
  w_upper <- rep(0,nMCMC)
  ### run L iterations for the chain lambda1, N1 and a1
  for (iter in 1:L){
    updateLambdaRes <- update_lambda(lambda1 - 1, partition1, ll_clusters1, p, V, dimV, a, N1)
    lambda1 <- updateLambdaRes$lambda + 1
    partition1 <- updateLambdaRes$clustering
    ll_clusters1 <- updateLambdaRes$clusterloglikelihoods
    ksize1 <- length(unique(lambda1))
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize1, g, n)
    N1 <- Ns[1]
    ### store 
    k1_chain[iter] <- ksize1
    N1_chain[iter] <- N1
    lambda1_chain[, iter] <- lambda1
  }
  ### try to couple the chains with lag L
  iter <- L + 1
  meet <- FALSE
  for (iter in (L + 1):nMCMC){
    ## update lambda 
    for (j in 1:n){
      lambdajs <- coupleLambda(j, lambda1, lambda2, p1 = p, p2 = p, a1 = a1, a2 = a2, N1 = N1, N2 = N2)
      # if(lambdajs[1] == lambdajs[2]){
      #   print(cat("for j = ", j ," coupled lambda-j = ", lambdajs[1]))
      # }
      lambda1[j] <- lambdajs[1]
      lambda2[j] <- lambdajs[2]
    }
    # ## relabel lambd1 and lambda2 and change a1 a2 accordingly
    relabel1 <- relabel(lambda1)
    lambda1 <- relabel1$lambda
    a1 <- a1[relabel1$iis,]
    relabel2 <- relabel(lambda2)
    lambda2 <- relabel2$lambda
    a2 <- a2[relabel2$iis,]
    lambda1_chain[, iter] <- lambda1
    lambda2_chain[, iter - L] <- lambda2 
    partition1 <- init_clustering(lambda1 - 1)
    partition2 <- init_clustering(lambda2 - 1)
    ksize1 <- length(unique(lambda1))
    ksize2 <- length(unique(lambda2))
    k1_chain[iter] <- ksize1
    k2_chain[iter - L] <- ksize2
    ## update alpha 
    ## update theta
    couple_theta <- coupleTheta(p1, partition1, p2, partition2)
    p1 <- couple_theta$p1
    p2 <- couple_theta$p2
    if(any(p1 == 0) | any(p2 == 0)){
      print(paste('couple_theta gives 0 at iteration',iter,'\n',sep = " "))
    }
    p1_chain[,iter] <- p1
    p2_chain[,iter - L] <- p2
    ## update N 
    Ns <- coupleN(N1 = N1, N2 = N2, k1 = ksize1, k2 = ksize2, g, n)
    N1 <- Ns[1]
    N2 <- Ns[2]
    N1_chain[iter] <- N1
    N2_chain[iter - L] <- N2
    percentageCoupled[iter] <- ( sum(lambda1 == lambda2) + sum(p1 == p2)) / (n + length(p1))
    if (iter %% (nMCMC / 50) == 0 ){
      print(paste('iteration', iter, percentageCoupled[iter] * 100, "percent coupled",sep = ' '))
    }
    if (N1 == N2 & all(lambda1 == lambda2) & all(p1 == p2) & !meet){
      meetTime <- iter
      meet <- TRUE
    }
    if (!meet & earlyStop){
      w_upper[iter - L] <- abs(N1_chain[iter] - N2_chain[iter - L]) +
        sum(abs(lambda1_chain[, iter] - lambda2_chain[, iter - L])) + 
        sum(abs(p1 - p2))
    }
    if (earlyStop & meet){
      break
    }
  }
  return(list(N1_chain = N1_chain, N2_chain = N2_chain,
              k1_chain = k1_chain, k2_chain = k2_chain,
              p1_chain = p1_chain, p2_chain = p2_chain, 
              meetTime = meetTime, percentageCoupled = percentageCoupled,
              lambda1_chain = lambda1_chain, lambda2_chain = lambda2_chain,
              L = L, w_upper = w_upper))
}

