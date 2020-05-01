# Coupling beta_0
# two functions copied from Pierre's github

get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)){
      return(list(xy = c(x,x), identical = TRUE))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq(1)
        reject <- (dq(y) + log(runif(1)) < dp(y))
      }
      return(list(xy = c(x,y), identical = FALSE))
    }
  }
}

rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
                        function(x) dnorm(x, mu1, sigma1, log = TRUE),
                        function(n) rnorm(n, mu2, sigma2),
                        function(x) dnorm(x, mu2, sigma2, log = TRUE))
  return(f())
}
#------------------------CoupleBeta0 (centred parametrisation)
#
# beta_prime_list_1, beta_prime_list_2: lists of length p containing beta_prime_j
# m_0, s_0_sq, s_sq: constants
#

coupleBeta0 <- function(beta_prime_list_1, beta_prime_list_2, mu_0, s_0_sq, s_sq){
  
  result <- lapply(1:length(beta_prime_list_1), function(l){
    k1 <- length(beta_prime_list_1[[l]])
    k2 <- length(beta_prime_list_2[[l]])
    
    if(k1==k2){
      # using the formula for the posterior in the Normal-Normal model
      common_var <- 1/(1/s_0_sq + k1/s_sq)
      mu1 <- common_var*(mu_0/s_0_sq + sum(beta_prime_list_1[[l]])/s_sq)
      mu2 <- common_var*(mu_0/s_0_sq + sum(beta_prime_list_2[[l]])/s_sq)
      return(rnorm_reflectionmxax(mu1, mu2, sqrt(common_var)))
      
    }else{
      var1 <- 1/(1/s_0_sq + k1/s_sq)
      var2 <- 1/(1/s_0_sq + k2/s_sq)
      mu1 <- var_1*(mu_0/s_0_sq + sum(beta_prime_list_1[[l]])/s_sq)
      mu2 <- var_2*(mu_0/s_0_sq + sum(beta_prime_list_2[[l]])/s_sq)
      return(rnorm_max_coupling(mu1, mu2, sqrt(var1), sqrt(var2)))
    }
  })
  # beta_0 for the first chain
  beta_0_1 <- unlist(lapply(result, function(x) x$xy[1]))
  # beta_0 for the second chain
  beta_0_2 <- unlist(lapply(result, function(x) x$xy[2]))
  beta_0_identical <- all(unlist(lapply(result, function(x) x$identical)))
  
  return(list(beta_0_1 = beta_0_1,
              beta_0_2 = beta_0_2,
              beta_0_identical = beta_0_identical))
  
}

#------------------------CoupleBeta0NonCentred (non-centred parametrisation)
# we are now working with differences between beta prime and beta_0 as suggested in the paper
#
# beta_prime_list_diffs_1, beta_prime_list_diffs_2: lists of length p containing beta_prime_j,l - beta_0,l
# previous_beta_0_1, previous_beta_0_2: previous values of beta_0_1, beta_0_2
# m_0, s_0_sq, s_sq: constants

coupleBeta0NonCentred <- function(beta_prime_list_diffs_1, beta_prime_list_diffs_2,
                                  previous_beta_0_1, previous_beta_0_2,
                                  mu_0, s_0_sq, s_sq){
  
  result <- lapply(1:length(beta_prime_list_diffs_1), function(l){
    k1 <- length(beta_prime_list_1[[l]])
    k2 <- length(beta_prime_list_2[[l]])
    
    if(k1==k2){
      # calculating beta_prime as previous_beta_0 + new differences
      beta_prime_1 <- previous_beta_0_1[l] + beta_prime_list_diffs_1[[l]]
      beta_prime_2 <- previous_beta_0_2[l] + beta_prime_list_diffs_2[[l]]
      
      # using the formula for the posterior in the Normal-Normal model
      common_var <- 1/(1/s_0_sq + k1/s_sq)
      mu1 <- common_var*(mu_0/s_0_sq + sum(beta_prime_1)/s_sq)
      mu2 <- common_var*(mu_0/s_0_sq + sum(beta_prime_2)/s_sq)
      return(rnorm_reflectionmxax(mu1, mu2, sqrt(common_var)))
      
    }else{
      # calculating beta_prime as previous_beta_0 + new differences
      beta_prime_1 <- previous_beta_0_1[l] + beta_prime_list_diffs_1[[l]]
      beta_prime_2 <- previous_beta_0_2[l] + beta_prime_list_diffs_2[[l]]
      
      # using the formula for the posterior in the Normal-Normal model
      # this time the variances may be different as they depend on the number of observations
      var1 <- 1/(1/s_0_sq + k1/s_sq)
      var2 <- 1/(1/s_0_sq + k2/s_sq)
      mu1 <- var_1*(mu_0/s_0_sq + sum(beta_prime_1)/s_sq)
      mu2 <- var_2*(mu_0/s_0_sq + sum(beta_prime_2)/s_sq)
      return(rnorm_max_coupling(mu1, mu2, sqrt(var1), sqrt(var2)))
    }
  })
  # beta_0 for the first chain
  beta_0_1 <- unlist(lapply(result, function(x) x$xy[1]))
  # beta_0 for the second chain
  beta_0_2 <- unlist(lapply(result, function(x) x$xy[2]))
  beta_0_identical <- all(unlist(lapply(result, function(x) x$identical)))
  
  return(list(beta_0_1 = beta_0_1,
              beta_0_2 = beta_0_2,
              beta_0_identical = beta_0_identical))
  
}


# analogous function as above but for a single kernel

SingleBeta0NonCentred <- function(beta_prime_list_diff,
                                  previous_beta_0,
                                  mu_0, s_0_sq, s_sq){
  
  result <- lapply(1:length(beta_prime_list_diff), function(l){
    
      # calculating beta_prime as previous_beta_0 + new differences
      beta_prime <- previous_beta_0[l] + beta_prime_list_diff[[l]]
      
      # using the formula for the posterior in the Normal-Normal model
      common_var <- 1/(1/s_0_sq + length(beta_prime_list_diff[[l]])/s_sq)
      mu <- common_var*(mu_0/s_0_sq + sum(beta_prime)/s_sq)
      return(rnorm(1,mu, sqrt(common_var)))
      
  })
  
  
  return(list(beta_0 = unlist(result)))
  
}

#---------------testing
#
# coupleBeta0(beta_prime_list_1 =list(rnorm(2), rnorm(3)), 
#             beta_prime_list_2 = list(rnorm(2), rnorm(2)),
#             mu_0 = 1, 
#             s_0_sq = 1, s_sq = 1)
# 

# SingleBeta0NonCentred(beta_prime_list_diffs = list(rnorm(10), rnorm(10)),
#                                   previous_beta = rnorm(2),
#                                   mu_0=1, s_0_sq=1, s_sq=1)
#   
