library(Rcpp)
source('couple_normal.R')
sourceCpp('single_kernel_alpha.R')


# ------------ update of beta0 (coupled kernel)
# recall we are working with differences between beta prime and beta_0 as suggested in the paper
#
# beta_diff_1, beta_diff_2 are matrices with n rows and p columns
# beta_0_1, beta_0_2 are a p-element vectors
# mu_0, s_0_sq, s_sq: hyperparameters

coupled_beta_0_update <- function(beta_diff_1,
                                  beta_diff_2,
                                  beta_0_1,
                                  beta_0_2,
                                  mu_0, 
                                  s_0_sq,
                                  s_sq,
                                  n){
  # using the formula for the posterior variance in the Normal-Normal model
  common_var <- 1/(1/s_0_sq + n/s_sq)
  
  result <- lapply(1:ncol(beta_diff), function(l){
    # calculating beta_prime as previous_beta_0 + new differences
    beta_prime_1 <- beta_0_1[l] + beta_diff_1[,l]
    beta_prime_2 <- beta_0_2[l] + beta_diff_2[,l]
    # posterior mean in the Normal-Normal model
    mu1 <- common_var*(mu_0/s_0_sq + sum(beta_prime_1)/s_sq)
    mu2 <- common_var*(mu_0/s_0_sq + sum(beta_prime_2)/s_sq)
    return(rnorm_reflectionmxax(mu1, mu2, sqrt(common_var)))
    
  })
  beta_0_1 <- unlist(lapply(result, function(x) x$xy[1]))
  beta_0_2 <- unlist(lapply(result, function(x) x$xy[2]))
  identical <- unlist(lapply(result, function(x) x$identical))
  return(list(beta_0_1 = beta_0_1, beta_0_2 = beta_0_2, identical = identical))
}


# ------------  function for performing an update of beta diff
# beta_diff_j_l is a value of beta_diff for cluster j and field l
# beta_0_l is a value of beta0 for field l
# l, icluster, clustering, theta_l, V are used as arguments in compute_loglikelihood_one_cluster_one_field_cpp
# partition_ll is an n x p matrix with current log-likelihood values for a given cluster and field
# s_sq is a hyperparameter
# proposal_sd is the std of the Gaussian proposal in the Metropolis-within-Gibbs update
# returns a list with the new beta_diff and the corresponding alpha for a given cluster and field

coupled_beta_diff_update <- function(beta_diff_j_l_1,
                                     beta_diff_j_l_2,
                                     beta_0_l_1,
                                     beta_0_l_2,
                                     l,
                                     icluster, 
                                     clustering_1,
                                     clustering_2,
                                     theta_l_1, 
                                     theta_l_2,
                                     V,
                                     s_sq,
                                     partition_ll_1,
                                     partition_ll_2,
                                     proposal_sd){
  # NA in this cluster means the cluster is empty so we are drawing from the prior
  if(is.na(partition_ll_1[icluster, l]) && is.na(partition_ll_2[icluster, l])){
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

# beta_diff_1, beta_diff_2 are an n x p matrices with value of beta diff for a given cluster (row) and field (column)
# alpha_1, alpha_2 are an n x p matrices with value of alpha for a given cluster (row) and field (column)
# beta_0_1, beta_0_2 are p-element vectors
# clustering_1, clustering_2, V are used as arguments in compute_loglikelihood_one_cluster_one_field_cpp
# theta_list_1, theta_list_2 are lists of vectors theta_l for l in 1:p used in compute_loglikelihood_one_cluster_one_field_cpp
# partition_ll_1, partition_ll_2 is an n x p matrix with current log-likelihood values for a given cluster and field
# mu_0, s_0_sq and s_sq are hyperparameters
# proposal_sd is the std of the Gaussian proposal in the Metropolis-within-Gibbs update
# p is the number of columns of V (variables)
# n is the number of rows of V (observations)
# returns a list with the the updated matrix beta_diff_1, beta_diff_2 and the corresponding matrix alpha_1, alpha_2 and upddated beta0_1, beta0_2

coupled_full_alpha_update <- function(beta_diff_1,
                                     beta_diff_2,
                                     alpha_1,
                                     alpha_2,
                                     beta_0_1,
                                     beta_0_2,
                                     clustering_1,
                                     clustering_2,
                                     theta_list_1,
                                     theta_list_2,
                                     V,
                                     partition_ll_1,
                                     partition_ll_2,
                                     mu_0, 
                                     s_0_sq,
                                     s_sq,
                                     proposal_sd,
                                     p,
                                     n){
  # updating beta_0 
  beta_0_update <- coupled_beta_0_update(beta_diff_1, beta_diff_2,
                                         beta_0_1, beta_0_2,
                                         mu_0, s_0_sq, s_sq, n)
  beta_0_1 <- beta_0_update$beta_0_1
  beta_0_2 <- beta_0_update$beta_0_2
  
  # matrices for storing new values
  new_beta_diff_1 <- matrix(NA, ncol = p, nrow = n)
  new_beta_diff_2 <- matrix(NA, ncol = p, nrow = n)
  new_alpha_1 <- matrix(NA, ncol = p, nrow = n)
  new_alpha_2 <- matrix(NA, ncol = p, nrow = n)
  # matrix for storing whether corresponding values of beta_diff are the same for the two chains
  new_beta_diff_identical <- matrix(NA, ncol = p, nrow = n)
  
  for(l in 1:p){
    theta_l_1 <- theta_list_1[[l]]
    theta_l_2 <- theta_list_2[[l]]
    beta_0_l_1 <- beta_0_1[l]
    beta_0_l_2 <- beta_0_2[l]
    for(icluster in 1:n){
      # update of beta prime and simultaneoulsy alpha
      # change this!!!!!
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
      
      new_beta_diff_1[icluster, l] <- result_list$beta_diff_j_l_1
      new_alpha_1[icluster, l] <- result_list$alpha_j_l_1
      new_beta_diff_2[icluster, l] <- result_list$beta_diff_j_l_2
      new_alpha_2[icluster, l] <- result_list$alpha_j_l_2
      new_beta_diff_identical[icluster, l] <- result_list$beta_diff_identical
    }
  }
  return(list(beta_diff_1 = new_beta_diff_1, beta_diff_2 = new_beta_diff_2,
              new_beta_diff_identical = new_beta_diff_identical,
              alpha_1 = new_alpha_1, alpha_2 = new_alpha_2,
              beta0_1 = beta_0_1, beta0_2 = beta_0_2, 
              beta_0_identical =   beta_0_update$identical))
  
}
