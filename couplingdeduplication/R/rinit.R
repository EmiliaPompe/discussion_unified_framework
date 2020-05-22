#'@title rinit
#'@description  initialization of the chains...
#'@param n the number of rows in the data
#'@param fieldfrequencies frequencies of appearance of each element in each field (list of vectors), 
#'such as provided by \code{\link{get_RLdata500}}
#'@param hyper a list of hyperparameters including m0, s0, s
#' where m0 is the mean of beta0, s0 its standard deviation, s is the std dev of beta_diff
#'@export
rinit <- function(n, fieldfrequencies, hyper, V){
  ## The state of the Markov chain  is (eta, N, beta_0, beta_diff, theta)
  p <- length(fieldfrequencies)
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
  ## initialize theta with frequencies of categories, initial values
  theta <- fieldfrequencies
  logtheta <- lapply(theta, function(l) log(l))
  # beta_0 <- rnorm(p, hyper$m0, hyper$s0)
  ## initial beta
  # beta_diff <- matrix(NA, nrow = n, ncol = p)
  # for (field in 1:p){
  #   beta_diff[,field] <- rnorm(n, 0, hyper$s)
  # }
  # alpha <- beta_to_alpha(beta_diff, beta_0)
  ## over-write alpha, beta_0 and beta_diff 
  alpha <- matrix(0.01, nrow = n, ncol = p)
  ## recall beta = logit(alpha) = log(alpha/(1-alpha))
  beta_0 <- rep(log(0.01/0.99), p)
  beta_diff <- matrix(0., nrow = n, ncol = p)
  ## compute log-likelihood associated with each cluster
  partition_ll <- 
    couplingdeduplication:::compute_loglikelihood_all_clusters_all_fields_cpp(partition, theta, logtheta, V-1, alpha)
  return(list(n = n, p = p, eta = eta, partition = partition, N = N, beta_0 = beta_0, beta_diff = beta_diff, 
              alpha = alpha, theta = theta, logtheta = logtheta, partition_ll = partition_ll))
}
