library(Rcpp)
sourceCpp('compute_loglikelihood_one_cluster_one_field.cpp')
# ------------ auxiliary function converting beta diff and beta0 to alpha
#
# beta_diff is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of beta_prime_{j,l} - beta_0_l corresponding to variable l and cluster j
# beta_0 is a p-element vector
# returns an n x p matrix alpha
beta_to_alpha <- function(beta_diff, beta_0){
  alpha <- beta_diff 
  for(i in 1:nrow(beta_diff)){
    
    exp_beta <- exp(beta_0 + beta_diff[i,])
    alpha[i,] <- exp_beta/(exp_beta+1)
    
  }
  return(alpha)
}


# ------------ auxiliary function converting alpha and beta0 to beta diff
#
# alpha is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of alpha corresponding to variable l and cluster j
# beta_0 is a p-element vector
# returns an n x p matrix beta_diff

alpha_to_beta <- function(alpha, beta_0){
  beta_diff <- alpha
  for(i in 1:nrow(alpha)){
    
    beta_diff[i,] <- log(alpha[i,]/(1-alpha[i,])) - beta_0
    
  }
  return(beta_diff)
}

# ------------ update of beta0 (single kernel)
# recall we are working with differences between beta prime and beta_0 as suggested in the paper
#
# beta_diff is a matrix with n rows and p columns
# beta_0 is a p-element vector
# m_0, s_0_sq, s_sq: hyperparameters

single_beta_0_update <- function(beta_diff,
                                 beta_0,
                                 mu_0, 
                                 s_0_sq,
                                 s_sq){
  # using the formula for the posterior variance in the Normal-Normal model
  common_var <- 1/(1/s_0_sq + nrow(beta_diff)/s_sq)
  
  result <- lapply(1:ncol(beta_diff), function(l){
    # calculating beta_prime as previous_beta_0 + new differences
    beta_prime <- beta_0[l] + beta_diff[,l]
    # posterior mean in the Normal-Normal model
    mu <- common_var*(mu_0/s_0_sq + sum(beta_prime)/s_sq)
    return(rnorm(1,mu, sqrt(common_var)))
    
  })
  return(list(beta_0 = unlist(result)))
}

# ------------ auxiliary function for calculating density for a given cluster and field

# the main argument is beta_diff_j_l: a value of beta_diff for cluster j and field l
# for calculating the prior we need  beta_diff_j_l and hyperparameter sq
# for calculating the likelihod we use the cpp function compute_loglikelihood_one_cluster_one_field_cpp
# this function takes l, icluster, clustering, theta_l, V and alpha_j_l (called a) as arguments
# alpha_j_l is alpha corresponding to a given beta_diff
# returns a value of the density for a given cluster and field

log_target_density_beta_prime <- function(beta_diff_j_l,
                                          l,
                                          icluster, 
                                          clustering,
                                          theta_l, 
                                          V,
                                          alpha_j_l,
                                          s_sq){
  
  # we use beta_diff hence the mean is 0, not beta_0
  logPrior <- dnorm(beta_diff_j_l, 0, sqrt(s_sq), log = TRUE)
  logLik <- compute_loglikelihood_one_cluster_one_field_cpp(l-1, icluster-1, clustering,
                                                            theta_l, V-1, alpha_j_l)
  return(logLik + logPrior)
}

# ------------  function for performing an update of beta diff
# beta_diff_j_l is a value of beta_diff for cluster j and field l
# beta_0_l is a value of beta0 for field l
# l, icluster, clustering, theta_l, V are used as arguments in compute_loglikelihood_one_cluster_one_field_cpp
# partition_ll is an n x p matrix with current log-likelihood values for a given cluster and field
# s_sq is a hyperparameter
# proposal_sd is the std of the Gaussian proposal in the Metropolis-within-Gibbs update
# returns a list with the new beta_diff and the corresponding alpha for a given cluster and field

single_beta_diff_update <- function(beta_diff_j_l,
                                    beta_0_l,
                                    l,
                                    icluster, 
                                    clustering,
                                    theta_l, 
                                    V,
                                    s_sq,
                                    partition_ll,
                                    proposal_sd){
  # NA in this cluster means the cluster is empty so we are drawing from the prior
  if(is.na(partition_ll[icluster, l])){
    beta_diff_j_l <- rnorm(1, 0, sqrt(s_sq))
    exp_beta <- exp(beta_0_l +beta_diff_j_l )
    alpha_j_l <- exp_beta/(exp_beta+1)
    return(list(beta_diff_j_l = beta_diff_j_l, alpha_j_l = alpha_j_l))
  }
  
  current_state <- beta_diff_j_l
  proposed_state <- rnorm(1, current_state, proposal_sd)
  # calculating the corresponding value of alpha
  exp_beta <- exp(beta_0_l + proposed_state)
  proposed_alpha_j_l <- exp_beta/(exp_beta+1)
  # common random number for acceptance/rejection
  logu <- log(runif(1))
  
  # lok likelihood taken from partition_ll + the prior
  current_density <-  partition_ll[icluster, l] + dnorm(beta_diff_j_l, 0, sqrt(s_sq), log = TRUE)
  proposed_density <- log_target_density_beta_prime(proposed_state, 
                                                    l, icluster, 
                                                    clustering, theta_l,
                                                    V, proposed_alpha_j_l,
                                                    s_sq)
  
  accept <- (logu < (proposed_density - current_density))
  
  if(accept){
    current_state <- proposed_state 
    alpha_j_l <- proposed_alpha_j_l
  }else{
    exp_beta <- exp(beta_0_l + current_state)
    alpha_j_l <- exp_beta/(exp_beta+1)
  }
  
  return(list(beta_diff_j_l = current_state, alpha_j_l = alpha_j_l))
}

# ------------  function for performing the full update on beta0, alpha and (equivalently) beta_diff

# beta_diff is an n x p matrix with value of beta diff for a given cluster (row) and field (column)
# alpha is an n x p matrix with value of alpha for a given cluster (row) and field (column)
# l, icluster, clustering, V are used as arguments in compute_loglikelihood_one_cluster_one_field_cpp
# theta_list is a list of vectors theta_l for l in 1:p used in compute_loglikelihood_one_cluster_one_field_cpp
# partition_ll is an n x p matrix with current log-likelihood values for a given cluster and field
# mu_0, s_0_sq and s_sq are hyperparameters
# proposal_sd is the std of the Gaussian proposal in the Metropolis-within-Gibbs update
# p is the number of columns of V (variables)
# n is the number of rows of V (observations)
# returns a list with the the updated matrix beta_diff and the corresponding matrix alpha, and  udpdated beta0

single_full_alpha_update <- function(beta_diff,
                                     alpha,
                                     beta_0,
                                     clustering,
                                     theta_list, 
                                     V,
                                     partition_ll,
                                     mu_0, 
                                     s_0_sq,
                                     s_sq,
                                     proposal_sd,
                                     p,
                                     n){
  # updating beta_0 
  beta_0_update <- single_beta_0_update(beta_diff, beta_0, mu_0, s_0_sq, s_sq)
  beta_0 <- beta_0_update$beta_0
  
  # matrices for storing new values
  new_beta_diff <- matrix(NA, ncol = p, nrow = n)
  new_alpha <- matrix(NA, ncol = p, nrow = n)
  
  for(l in 1:p){
    theta_l <- theta_list[[l]]
    beta_0_l <- beta_0[l]
    for(icluster in 1:n){
      # update of beta prime and simultaneoulsy alpha
      result_list <- single_beta_diff_update(beta_diff_j_l = beta_diff[icluster, l],
                                             beta_0_l = beta_0_l,
                                             l = l,
                                             icluster = icluster, 
                                             clustering = clustering,
                                             theta_l = theta_l, 
                                             V = V,
                                             s_sq = s_sq,
                                             partition_ll = partition_ll,
                                             proposal_sd = proposal_sd)
      new_beta_diff[icluster, l] <- result_list$beta_diff_j_l
      new_alpha[icluster, l] <- result_list$alpha_j_l
    }
  }
  return(list(beta_diff = new_beta_diff, alpha = new_alpha, beta0 = beta_0))
  
}
