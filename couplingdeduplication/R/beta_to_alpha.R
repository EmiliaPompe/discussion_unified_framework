
# ------------ auxiliary function converting beta diff and beta_0 to alpha
#
# beta_diff is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of beta_prime_{j,l} - beta_0_l corresponding to variable l and cluster j
# beta_0 is a p-element vector
# returns an n x p matrix alpha
#'@export
beta_to_alpha <- function(beta_diff, beta_0){
  alpha <- beta_diff 
  for(i in 1:nrow(beta_diff)){
    
    exp_beta <- exp(beta_0 + beta_diff[i,])
    alpha[i,] <- exp_beta/(exp_beta+1)
    
  }
  return(alpha)
}


