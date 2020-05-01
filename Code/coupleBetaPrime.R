source('couple_normal.R')

# Coupling beta prime
# for a fixed l = 1,...,p and a fixed cluster j
# (this makes us work with univariate variables only)

# beta is the univariate variable of interest
# v_vector is the vector of values v for the fixed cluster and variable
# theta_vector gives the corresponding values of theta for the vector v above
# beta_0 gives the prior mean
# s_sq is the prior variance


logTargetDensityBetaPrime <- function(beta, theta_vector, v_vector, beta_0, s_sq){
  logPrior <- dnorm(beta, beta_0, sqrt(s_sq))
  # defining alpha as in the bottom of page 2 of our write-up
  alpha <- exp(beta)/(1+ exp(beta))

  # the recursion formula (6) from our write-up
  recursion_formula_vector <- rep(NA, length(v_vector))
  recursion_formula_vector[1] <- theta_vector[1]
  if(length(v_vector)>1){
    for(k in 2:length(v_vector)){
      
      aux_vector <- sapply(1:(k-1), function(h){
        aux <- alpha*theta_vector[h]
        if(v_vector[k]==v_vector[h]){
          aux <- aux + (1-alpha)
        }
        return(aux)
      })
      recursion_formula_vector[k] <- theta_vector[k]*(alpha*recursion_formula_vector[k-1] + (1-alpha)*prod(aux_vector))
    }
  }
  # log likelihood is the last element of the recursive formula
  logLik <- log(tail(recursion_formula_vector,1))
  logLik + logPrior
}

# current_beta_prime_1 and current_beta_prime_2 are current values of beta_prime for both chains
# current_beta_prime_ident is TRUE/FALSE defining if the above variables are equal
# proposal_sd is the proposal standard deviation used for proposing a new point
# the remaining variables are as above for both chains

CoupleBetaPrime <- function(current_beta_prime_1, current_beta_prime_2,
                            current_beta_prime_ident,
                            theta_vector_1, theta_vector_2, 
                            v_vector_1, v_vector_2, 
                            beta_0_1, beta_0_2,
                            s_sq, proposal_sd){
  logTargetDensity1 <- function(beta){
    logTargetDensityBetaPrime(beta, theta_vector_1, v_vector_1, beta_0_1, s_sq)
  }
  logTargetDensity2 <- function(beta){
    logTargetDensityBetaPrime(beta, theta_vector_2, v_vector_2, beta_0_2, s_sq)
  }
  
  # drawing proposals from maximal coupling
  max_coupling <- rnorm_reflectionmxax(current_beta_prime_1, current_beta_prime_2, proposal_sd)
  
  # common random number for acceptance/rejection
  logu <- log(runif(1))
  current_density_1 <- logTargetDensity1(current_beta_prime_1)
  current_density_2 <- logTargetDensity2(current_beta_prime_2)
  proposed_density_1 <- logTargetDensity1(max_coupling$xy[1])
  current_density_2 <- logTargetDensity2(max_coupling$xy[2])
  
  accept1 <- (logu < (proposed_density_1 - current_density_1))
  if(accept1){
    beta_prime_1 <- max_coupling$xy[1] 
  }else{
    beta_prime_1 <- max_coupling$xy[1]
  }
  accept2 <- (logu < (proposed_density_1 - current_density_1))

  # checking if the coupling has already happened
  ident <- max_coupling$identical && accept1 && accept2
  
  # in case they were previously coupled and now both moves are rejected
  if((!ident) && current_beta_prime_ident){
    ident <- (!accept1) && (!accept2)
  }
  
  return(list(beta_prime_1 = beta_prime_1,
             beta_prime_2 = beta_prime_2,
             beta_prime_identical = ident))
}

#------------------------logTargetDensityBetaPrimeNonCentred (non-centred parametrisation)
#
# we are now working with differences between beta prime and beta_0 denoted by beta_diff
# beta_diff is a difference beta_prime_{j,l} - beta_{0,l}
# v_vector is the vector of values v for the fixed cluster and variable
# theta_vector gives the corresponding values of theta for the vector v above
# beta_0 gives the prior mean
# s_sq is the prior variance

logTargetDensityBetaPrimeNonCentred <- function(beta_diff, theta_vector, v_vector, beta_0, s_sq){
  beta <- beta_0 + beta_diff
  logPrior <- dnorm(beta, beta_0, sqrt(s_sq))
  # defining alpha as in the bottom of page 2 of our write-up
  alpha <- exp(beta)/(1+ exp(beta))
  
  # the recursion formula (6) from our write-up
  recursion_formula_vector <- rep(NA, length(v_vector))
  recursion_formula_vector[1] <- theta_vector[1]
  if(length(v_vector)>1){
    for(k in 2:length(v_vector)){
      
      aux_vector <- sapply(1:(k-1), function(h){
        aux <- alpha*theta_vector[h]
        if(v_vector[k]==v_vector[h]){
          aux <- aux + (1-alpha)
        }
        return(aux)
      })
      recursion_formula_vector[k] <- theta_vector[k]*(alpha*recursion_formula_vector[k-1] + (1-alpha)*prod(aux_vector))
    }
  }
  # log likelihood is the last element of the recursive formula
  logLik <- log(tail(recursion_formula_vector,1))
  logLik + logPrior
}

#------------------------CoupleBetaPrimeNonCentred (non-centred parametrisation)
#
# current_beta_prime_diff_1 and current_beta_prime_diff_2 are current values of differences beta_prime_{j,l} - beta_{0,l} for both chains
# (hence current_beta_prime_diff_1 and current_beta_prime_diff_2 are vectors of length 1)
# current_beta_prime_diff_ident is TRUE/FALSE defining if the above variables are equal
# proposal_sd is the proposal standard deviation used for proposing a new point
# the remaining variables are as above for both chains
CoupleBetaPrimeNonCentred <- function(current_beta_prime_diff_1, current_beta_prime_diff_2,
                                      current_beta_prime_diff_ident,
                                      theta_vector_1, theta_vector_2, 
                                      v_vector_1, v_vector_2, 
                                      beta_0_1, beta_0_2,
                                      s_sq, proposal_sd){
  logTargetDensity1 <- function(beta_diff){
    logTargetDensityBetaPrimeNonCentred(beta_diff, theta_vector_1, v_vector_1, beta_0_1, s_sq)
  }
  logTargetDensity2 <- function(beta_diff){
    logTargetDensityBetaPrimeNonCentred(beta_diff, theta_vector_2, v_vector_2, beta_0_2, s_sq)
  }
  
  # drawing proposals from maximal coupling
  max_coupling <- rnorm_reflectionmxax(current_beta_prime_diff_1, current_beta_prime_diff_2, proposal_sd)
  
  # common random number for acceptance/rejection
  logu <- log(runif(1))
  current_density_1 <- logTargetDensity1(current_beta_prime_diff_1)
  current_density_2 <- logTargetDensity2(current_beta_prime_diff_2)
  proposed_density_1 <- logTargetDensity1(max_coupling$xy[1])
  current_density_2 <- logTargetDensity2(max_coupling$xy[2])
  
  accept1 <- (logu < (proposed_density_1 - current_density_1))
  if(accept1){
    beta_prime_1 <- max_coupling$xy[1] 
  }else{
    beta_prime_1 <- max_coupling$xy[1]
  }
  accept2 <- (logu < (proposed_density_1 - current_density_1))
  
  # checking if the coupling has already happened
  ident <- max_coupling$identical && accept1 && accept2
  
  # in case they were previously coupled and now both moves are rejected
  if((!ident) && current_beta_prime_diff_ident){
    ident <- (!accept1) && (!accept2)
  }
  
  return(list(beta_prime_diff_1 = beta_prime_1,
              beta_prime_diff_2 = beta_prime_2,
              beta_prime_identical = ident))
}


# analogous function as above but for a single kernel
SingleBetaPrimeNonCentred <- function(current_beta_prime_diff, 
                                      theta_vector, 
                                      v_vector, 
                                      beta_0,
                                      s_sq, proposal_sd){
  logTargetDensity <- function(beta_diff){
    logTargetDensityBetaPrimeNonCentred(beta_diff, theta_vector, v_vector, beta_0, s_sq)
  }
  
  current_state <- current_beta_prime_diff
  # common random number for acceptance/rejection
  logu <- log(runif(1))
  current_density <- logTargetDensity(current_state)
  proposed_state <- rnorm(1, current_state, proposal_sd)
  proposed_density <- logTargetDensity(proposed_state)
  
  
  accept <- (logu < (proposed_density - current_density))
  if(accept){
   current_state <- proposed_state  
  }
  
  return(list(beta_prime_diff = current_state))
}
#--------------------- testing
#
# SingleBetaPrimeNonCentred(current_beta_prime_diff = rnorm(2), 
#                                       theta_vector = runif(2), 
#                                       v_vector = c(3,4), 
#                                       beta_0 = 2.5,
#                                       s_sq = 0.1, proposal_sd = 0.1)
# 
