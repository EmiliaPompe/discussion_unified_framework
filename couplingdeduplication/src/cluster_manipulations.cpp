#include "cluster_manipulations.h"
using namespace Rcpp;
using namespace std;

void remove_label_from_partition(int ieta, int label,
                                 int & ksize,
                                 IntegerMatrix & clmembers,
                                 IntegerVector & clsize,
                                 NumericMatrix & clusterloglikelihoods,
                                 const List & theta,
                                 const List & logtheta,
                                 const IntegerMatrix & V,
                                 const NumericMatrix & alpha){
  int n = V.nrow();
  int p = V.ncol();
  // i.e. translate -1's towards the left
  bool shift = false;
  for (int icol = 0; icol < clsize[label]; icol++){
    if (clmembers(label,icol) == ieta){
      shift = true;
    }
    if (shift){
      // move entry of next column to current column
      if (icol + 1 < n){
        clmembers(label,icol) = clmembers(label,icol+1);
      } else {
        clmembers(label,icol) = -1;
      }
    } 
  }
  clsize[label] --;
  // if cluster is now empty
  if (clsize[label] == 0){
    // decrement cluster counter
    ksize --;
    // erase log-likelihoods associated with that cluster
    std::fill(clusterloglikelihoods.row(label).begin(), clusterloglikelihoods.row(label).end(), NA_REAL);
  } else {
    // recompute likelihood of cluster from scratch
    for (int field = 0; field < p; field ++){
      clusterloglikelihoods(label,field) = ll_cluster_field(clsize[label],
                            clmembers(label,_), theta[field], logtheta[field], V(_,field), alpha(label,field));
    }        
  }
}

void compute_proba_eta(NumericMatrix & uponetajoining_loglikelihood, 
                        NumericVector & logproba_eta, 
                        NumericVector & newblock_loglikelihood, NumericVector & theta_field, NumericVector & logtheta_field, 
                        int ieta, int p, int n, 
                        const IntegerVector & clsize,
                        const IntegerMatrix & clmembers,
                        const NumericMatrix & clusterloglikelihoods,
                        const List & theta, const List & logtheta,
                        const IntegerMatrix & V,
                        const NumericMatrix & alpha,
                        const int  N, const int ksize){
  // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
  // next we can compute the likelihood associated with a new block with that label in it
  for (int field = 0; field < p; field++){
    logtheta_field = logtheta[field];
    newblock_loglikelihood(field) = logtheta_field(V(ieta,field));
  }
  // vector to aggregate likelihood and prior probabilities for new label
  std::fill(logproba_eta.begin(), logproba_eta.end(), 0.);
  // next we loop over clusters
  for (int icluster = 0; icluster < n; icluster ++){
    // first compute likelihoods of joining
    if (clsize[icluster]==0){
      // this would be a new cluster
      uponetajoining_loglikelihood.row(icluster) = newblock_loglikelihood;
    } else {
      // this would be an existing cluster
      for (int field = 0; field < p; field++){
        theta_field = theta[field];
        logtheta_field = logtheta[field];
        // first part of recursion
        uponetajoining_loglikelihood(icluster, field) = log(alpha(icluster,field)) + logtheta_field(V(ieta,field)) + 
          clusterloglikelihoods(icluster,field);
        double logprod = 0.;
        // next, implement second part of recursion
        for (int clustermember = 0; clustermember < clsize[icluster]; clustermember ++){
          int qprime = clmembers(icluster, clustermember);
          logprod += log((1 - alpha(icluster,field)) * (V(ieta,field) == V(qprime,field)) + alpha(icluster,field) * theta_field(V(qprime,field)));
        }
        logprod += log(1 - alpha(icluster,field)) + logtheta_field(V(ieta,field));
        double max_logs = std::max(logprod, uponetajoining_loglikelihood(icluster, field));
        uponetajoining_loglikelihood(icluster, field) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood(icluster, field) - max_logs));
      }
    }
    // then aggregate priors and conditional likelihood
    if (clsize[icluster] == 0){
      // i.e. P(eta[ieta] = q) where q is the label of a new block
      logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster));
      logproba_eta[icluster] += (log(N-ksize) - log(n-ksize)); 
      // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster 
    } else {
      logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster) - clusterloglikelihoods.row(icluster));
      logproba_eta[icluster] += 0.;
    }
  }
}
