### this function relabels the lambda's such that partitions are fixed
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


# lambda <- c(1,4,3,4,2)
relabel <- function(lambda){
  ## returns new labels and the correspondence between old labels and new labels 
  n <- length(lambda)
  iis <- c(1:n)
  newlambda <- rep(NA, n)
  visited <- rep(FALSE, n)
  changed <- rep(FALSE, n)
  currentCnt <- 1
  for ( j in 1:n){
    if(!visited[j]){
      iis[lambda[j]] <- currentCnt
      changed[lambda[j]] <- TRUE
      clMembers <- which(lambda == lambda[j])
      visited[clMembers] <- TRUE
      newlambda[clMembers] <- currentCnt
      currentCnt <- currentCnt + 1
    }
  }
  iis[!changed] <- c(currentCnt : n)
  return(list(lambda = newlambda, iis = iis))
}
