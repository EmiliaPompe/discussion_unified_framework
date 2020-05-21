# function converting alpha and beta_0 to beta diff
# alpha is a matrix with n rows and p columns; each row denotes a cluster and each column denotes a variable
# if a cluster is empty we keep NA throughout the row; otherwise we keep the value of alpha corresponding to variable l and cluster j
# beta_0 is a p-element vector
# returns an n x p matrix beta_diff
#'@export
alpha_to_beta <- function(alpha, beta_0){
  beta_diff <- alpha
  for(i in 1:nrow(alpha)){
    beta_diff[i,] <- log(alpha[i,]/(1-alpha[i,])) - beta_0
  }
  return(beta_diff)
}
