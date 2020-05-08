### this function relabels the eta's such that partitions are fixed
### but labels are sorted according to first appearance
### for example:
# One simple re-labeling would be to re-assign the labels by increasing order 
# with respect to the record indices in the data set. After each sweep of eta-updates, 
# for example, if we have 
# eta = {1,4,3,4,2}
# and associated u'_z = {1,4,3,2}, Z = {1|24|3|5}
# we can relabel eta into
# eta = {1,2,3,2,4} (increasing order of first appearance for each value of eta)
# and associated u'_z = {1,2,3,4} and Z = {1|24|3|5} as before.

relabel <- function(eta, clustering){
  ## returns new labels and the correspondence between old labels and new labels 
  n <- length(eta)
  old_to_new <- c(1:n) ## maps old cluster index to new cluster index
  new_eta <- rep(NA, n)
  visited <- rep(FALSE, n)
  changed <- rep(FALSE, n)
  currentCnt <- 1
  for (j in 1 : n){
    if(!visited[j]){
      old_to_new[eta[j]] <- currentCnt
      changed[eta[j]] <- TRUE
      clMembers <- which(eta == eta[j])
      visited[clMembers] <- TRUE
      new_eta[clMembers] <- currentCnt
      currentCnt <- currentCnt + 1
    }
  }
  old_to_new[!changed] <- c(currentCnt : n)
  return(list(eta = new_eta, old_to_new = old_to_new))
}

## function to perform re-labeling
## this one goes through each cluster and gives it the label
## of smallest index of member of cluster
relabel2 <- function(eta, clustering){
  new_eta <- rep(0, length(eta))
  ## permutation is a vector with i-th element equal to the new label given to what 
  ## was previously the i-th label 
  permutation <- rep(0, length(eta))
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
  permutation[permutation==0] <- setdiff(1:length(eta), sort(permutation[permutation!=0]))
  return(list(eta = new_eta, old_to_new = order(permutation)))
  # return(list(eta = new_eta, clustering = list(clsize = clustering$clsize[order(permutation)], 
  #                                              ksize = clustering$ksize,
  #                                              clmembers = clustering$clmembers[order(permutation),])))
}

