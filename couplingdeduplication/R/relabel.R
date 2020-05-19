## function to perform re-labeling
## goes through each cluster and gives it the label
## of smallest index among members of the cluster
#'@export
relabel <- function(eta, clustering){
  n <- length(eta)
  new_eta <- rep(0, n)
  ## permutation is a vector with i-th element equal to the new label given to what 
  ## was previously the i-th label 
  permutation <- rep(0, n)
  ## for each cluster...
  for (icluster in 1:length(clustering$clsize)){
    ## if cluster is not empty
    if (clustering$clsize[icluster]>0){
      ## get indices of cluster members
      cluster_members <- clustering$clmembers[icluster,1:clustering$clsize[icluster]]
      ## give as label smallest of these indices
      cluster_label <- min(cluster_members)
      new_eta[cluster_members+1] <- cluster_label + 1
      permutation[icluster] <- cluster_label + 1
    }
  }
  ## assign arbitrary labels to empty clusters
  permutation[permutation==0] <- setdiff(1:n, sort(permutation[permutation!=0]))
  return(list(eta = new_eta, old_to_new = order(permutation)))
}

