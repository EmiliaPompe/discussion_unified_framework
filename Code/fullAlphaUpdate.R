source('coupleBeta0.R')
source('coupleBetaPrime.R')

# ------------Auxiliary function converting alpha notation to beta prime
#
# alpha_matrix is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of alpha corresponding to variable l and cluster j
# beta_0 is a p-element vector

AlphaToBetaPrime <- function(alpha_matrix, beta_0){
  beta_prime_diff <- alpha_matrix
  for(i in 1:nrow(alpha_matrix)){
    if(!is.na(alpha_matrix[i,1])){
      beta_prime_diff[i,] <- log(alpha_matrix[i,]/(1-alpha_matrix[i,])) - beta_0
    }
  }
  return(beta_prime_diff)
}

# ------------Auxiliary function converting beta prime notation to alpha
#
# beta_prime_diff_matrix is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of beta_prime_{j,l} - beta_0_l corresponding to variable l and cluster j
# beta_0 is a p-element vector
BetaPrimeToAlpha <- function(beta_prime_diff_matrix, beta_0){
  alpha_matrix <- beta_prime_diff_matrix 
  for(i in 1:nrow(beta_prime_diff_matrix)){
    if(!is.na(beta_prime_diff_matrix[i,1])){
      exp_beta <- exp(beta_0 + beta_prime_diff_matrix[i,])
      alpha_matrix[i,] <- exp_beta/(exp_beta+1)
    }
  }
  return(alpha_matrix)
}

#-------------- Function for performing a full update of alpha for a single kernel
# alpha_matrix and beta_0 defined as above
# proposal_sd, mu_0, s_0_sq, s_sq defined as in CoupleBetaPrime.R
# theta_list is a p-element list of vectors denoting the probability of values of v for l-th variable
# V_matrix is the n x p matrix of observations where each row corresponds to a record and each column to a variable
# eta is the vector of cluster indices as used earlier

singleAlphaUpdate <- function(alpha_matrix,
                              beta_0,
                              theta_list,
                              V_matrix,
                              eta,
                              proposal_sd,
                              mu_0, s_0_sq, s_sq){
  # converting the matrix to beta prime values
  beta_prime_diff_matrix <- AlphaToBetaPrime(alpha_matrix, beta_0)
  # p is the number of variables
  p <- ncol(V_matrix) 
  # what clusters currently are non-empty
  current_clusters <- which(!is.na(beta_prime_diff_matrix[,1]))
  beta_diff_list <- lapply(1:p, function(l) beta_prime_diff_matrix[current_clusters,l])
  # performing the update of beta0
  update_beta_0 <- SingleBeta0NonCentred(beta_diff_list,
                                         beta_0,
                                         mu_0, s_0_sq, s_sq)
  current_beta_0 <- update_beta_0$beta_0
  # the list below stores indices of observations belonging to the clusters we currently have
  cluster_indices <- lapply(current_clusters, function(cluster_nr) which(eta == cluster_nr))
  
  result <- lapply(1:p, function(l){
    beta_prime_diff_l <- lapply(1:length(current_clusters), function(k){
      # values of v for a given variable and cluster
      v_vector <- V_matrix[cluster_indices[[k]],l]
      # values of theta for corresponding to v_vector and variable l
      theta_vector<- theta_list[[l]][v_vector]
      # updating beta_prime diff
      beta_prime_diff_update <- SingleBetaPrimeNonCentred(beta_diff_list[[l]][k], 
                                                          theta_vector, 
                                                          v_vector, 
                                                          current_beta_0[l],
                                                          s_sq, proposal_sd)
      return(beta_prime_diff_update$beta_prime_diff)
    })
    return(unlist(beta_prime_diff_l))
  })
  # updating beta prime diff  
  beta_prime_diff_matrix[current_clusters, 1:p] <- matrix(unlist(result), ncol = p)
  # we convert it back to alpha 
  return(BetaPrimeToAlpha(beta_prime_diff_matrix, current_beta_0))
}




#-------------- Function for performing a full update of alpha for a coupled kernel
# alpha_matrix_1, alpha_matrix_2: matrices with n rows and p columns; each row denotes a cluster and each column denotes a variable
# theta_list is a p-element list of vectors denoting the probability of values of v for l-th variable
# beta_0_1, beta_0_2 are p-element vectors
# V_matrix is the n x p matrix of observations where each row corresponds to a record and each column to a variable
# proposal_sd, mu_0, s_0_sq, s_sq defined as in CoupleBetaPrime.R


# eta is the vector of cluster indices as used earlier

coupledAlphaUpdate <- function(alpha_matrix_1,
                               alpha_matrix_2,
                               alpha_identical,
                               beta_0_1,
                               beta_0_2,
                               theta_list_1,
                               theta_list_2,
                               V_matrix,
                               eta_1,
                               eta_2,
                               proposal_sd,
                               mu_0, s_0_sq, s_sq){
  # p is the number of variables
  p <- ncol(V_matrix) 
  
  # converting the matrix to beta prime values
  beta_prime_diff_matrix_1 <- AlphaToBetaPrime(alpha_matrix_1, beta_0_1)
  beta_prime_diff_matrix_2 <- AlphaToBetaPrime(alpha_matrix_2, beta_0_2)
  # what clusters currently are non-empty
  current_clusters_1 <- which(!is.na(beta_prime_diff_matrix_1[,1]))
  current_clusters_2 <- which(!is.na(beta_prime_diff_matrix_2[,1]))
  
  # which clusters are only in chain 1, in chain 2 and in both chains
  clusters_both_chains <- intersect(current_clusters_1, current_clusters_2)
  clusters_only_1 <- setdiff(current_clusters_1, current_clusters_2)
  clusters_only_2 <- setdiff(current_clusters_2, current_clusters_1)
  
  
  full_cluster_indices_1 <- lapply(1:nrow(alpha_matrix_1), function(cluster_nr) which(eta_1 == cluster_nr))
  full_cluster_indices_2 <- lapply(1:nrow(alpha_matrix_2), function(cluster_nr) which(eta_2 == cluster_nr))
  
  beta_diff_list_1 <- lapply(1:p, function(l) beta_prime_diff_matrix_1[current_clusters_1,l])
  beta_diff_list_2 <- lapply(1:p, function(l) beta_prime_diff_matrix_2[current_clusters_2,l])
  
  # coupled update of beta_0
  update_beta_0 <- coupleBeta0(beta_diff_list_1, beta_diff_list_2, mu_0, s_0_sq, s_sq)
  current_beta_0_1 <- update_beta_0$beta_0_1
  current_beta_0_2 <- update_beta_0$beta_0_2
  beta_0_identical <- update_beta_0$beta_0_identical
  
  
  # updating beta prime diff  
  current_beta_prime_diff_matrix_1 <- matrix(NA, ncol = p, nrow = nrow(alpha_matrix_1))
  current_beta_prime_diff_matrix_2 <- matrix(NA, ncol = p, nrow = nrow(alpha_matrix_1))
  current_alpha_identical <- matrix(NA, ncol = p, nrow = nrow(alpha_matrix_1))
  
  for(l in 1:p){
    # has beta_0 for variable l coupled?
    beta_0_l_ident <- beta_0_identical[l]
    
    for(j in 1:nrow(alpha_matrix_1)){
      # we look separately at cases when cluster j is empty in only one of the chains, non-empty in both or empty in both
      # in case when only one cluster in non-empty, we use a non-coupled kernel for the update
      
      if(j %in% clusters_only_1){
        # values of v for a given variable and cluster
        v_vector <- V_matrix[ full_cluster_indices_1[[j]],l]
        # values of theta for corresponding to v_vector and variable l
        theta_vector<- theta_list_1[[l]][v_vector]
        # updating beta_prime diff
        beta_prime_diff_update <- SingleBetaPrimeNonCentred(beta_prime_diff_matrix_1[j,l], 
                                                            theta_vector, 
                                                            v_vector, 
                                                            current_beta_0_1[l],
                                                            s_sq, proposal_sd)
        current_list <- list(beta_prime_diff_update$beta_prime_diff, NA, FALSE)
      } else if(j %in% clusters_only_2){
        # values of v for a given variable and cluster
        v_vector <- V_matrix[ full_cluster_indices_2[[j]],l]
        # values of theta for corresponding to v_vector and variable l
        theta_vector<- theta_list_2[[l]][v_vector]
        # updating beta_prime diff
        beta_prime_diff_update <- SingleBetaPrimeNonCentred(beta_prime_diff_matrix_2[j,l], 
                                                            theta_vector, 
                                                            v_vector, 
                                                            current_beta_0_2[l],
                                                            s_sq, proposal_sd)
        
        return(list(NA, beta_prime_diff_update$beta_prime_diff, FALSE))
        
      } else if(j %in% clusters_both_chains){
        
        # values of v for a given variable and cluster for both chains
        v_vector_1 <- V_matrix[ full_cluster_indices_1[[j]],l]
        v_vector_2 <- V_matrix[ full_cluster_indices_2[[j]],l]
        # values of theta for corresponding to v_vector and variable l for both chains
        theta_vector_1<- theta_list_1[[l]][v_vector_1]
        theta_vector_2<- theta_list_2[[l]][v_vector_2]
        
        beta_prime_diff_update <- CoupleBetaPrimeNonCentred(beta_prime_diff_matrix_1[j,l],
                                                            beta_prime_diff_matrix_2[j,l],
                                                            alpha_identical[j,l],
                                                            theta_vector_1,
                                                            theta_vector_2, 
                                                            v_vector_1,
                                                            v_vector_2, 
                                                            current_beta_0_1[l],
                                                            current_beta_0_2[l],
                                                            s_sq, proposal_sd)
        current_list <- list(beta_prime_diff_update$beta_prime_diff_1, 
                             beta_prime_diff_update$beta_prime_diff_2,
                             beta_prime_diff_update$beta_prime_identical)
        
        
      } else{
        # we assume that if cluster j was empty for both chains than beta prime's have coupled; see the explanation for alpha below
       current_list <- list(NA, NA, TRUE)
      }
      current_beta_prime_diff_matrix_1[j,l] <- current_list[[1]]
      current_beta_prime_diff_matrix_2[j,l] <- current_list[[2]]
      # the values of alpha are identical in both chains iff both beta_0 and beta prime's have coupled
      # empty clusters can be thought of as driven exclusively by beta_0 hence we have current_list[[3]] = TRUE in this scenario
      current_alpha_identical[j,l] <- current_list[[3]] && beta_0_l_ident
    }
  }
  
  
  
}


# ----------- testing
p <- 3
alpha_matrix <- matrix(c(runif(6), NA, NA, NA), ncol =p, byrow = T)
beta_0 <- rnorm(3, 0)
beta_prime_diff_matrix <- AlphaToBetaPrime(alpha_matrix, beta_0)
alpha_matrix_2 <- BetaPrimeToAlpha(beta_prime_diff_matrix, beta_0)
all.equal(alpha_matrix, alpha_matrix_2)

eta <- c(rep(1,4), rep(2,5), 1)
V_matrix <- matrix(sample(1:3, 30,  replace = TRUE), ncol =3)
theta_list <- list(c(0.4, 0.4, 0.2), c(0.4, 0.5, 0.1), c(0.2, 0.4, 0.4))

singleAlphaUpdate(alpha_matrix,
                  beta_0,
                  theta_list,
                  V_matrix,
                  eta,
                  0.1,
                  1, 1, 1)



