## Initialization of the chains
## n the number of rows in the data
## fieldfrequencies frequencies of appearance of each element in each field (list of vectors)
## hyper a list of hyperparameters including m0, s0, s
#'@export
rinit <- function(n, fieldfrequencies, hyper){
  ## The state of the Markov chain  is (eta, N, beta_0, beta_diff, theta)
  ## initial eta 
  eta <- sample(x = 1:n, size = n, replace = T)
  ## initial N
  N <- 2500
  ## initial clustering based on eta
  partition <- couplingdeduplication:::init_clustering_cpp(eta-1)
  ## relabel clustering
  relabel_result <- relabel(eta, partition)
  eta <- relabel_result$eta
  partition$clsize <- partition$clsize[relabel_result$old_to_new] 
  partition$clmembers <- partition$clmembers[relabel_result$old_to_new,]
  ## initial b0
  p <- length(fieldfrequencies)
  beta_0 <- rnorm(p, hyper$m0, hyper$s0)
  ## initial beta
  beta_diff <- matrix(NA, nrow = n, ncol = p)
  for (field in 1:p){
    beta_diff[,field] <- rnorm(n, 0, hyper$s)
  }
  ## initialize theta with frequencies of categories, initial values
  theta <- fieldfrequencies
  ## these are deterministic transformations, to save time later (?)
  alpha <- beta_to_alpha(beta_diff, beta_0)
  logtheta <- lapply(theta, function(l) log(l))
  
  return(list(n = n, p = p, eta = eta, partition = partition, N = N, beta_0 = beta_0, beta_diff = beta_diff, 
              alpha = alpha, theta = theta, logtheta = logtheta))
}
