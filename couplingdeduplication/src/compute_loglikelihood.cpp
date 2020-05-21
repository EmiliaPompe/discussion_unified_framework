#include <Rcpp.h>
#include "compute_loglikelihood.h"
#include "ll_cluster_field.h"
using namespace Rcpp;

// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times H = ncol(V)
// if a cluster is empty, the corresponding row is filled with zeros
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to 
// values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// a is alpha_ prime[icluster,l] in the paper, and it is alpha[icluster,l] of
// the R code

// [[Rcpp::export]]
double compute_loglikelihood_one_cluster_one_field_cpp(int l ,
                                                       int icluster, 
                                                       const List & clustering,
                                                       const NumericVector & theta_l,
                                                       const NumericVector & logtheta_l,
                                                       const IntegerMatrix & V,
                                                       const double a) {
  const IntegerVector & clsize = clustering["clsize"];
  const IntegerMatrix & clmembers = clustering["clmembers"];
  double log_likelihood = NA_REAL ; 
  log_likelihood = ll_cluster_field(clsize[icluster], 
                                    clmembers(icluster,_),
                                    theta_l, 
                                    logtheta_l,
                                    V(_,l), a);
  return log_likelihood;
}


// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with NA
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// alpha is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericMatrix compute_loglikelihood_all_clusters_all_fields_cpp(const List & clustering,
                                                                const List & theta,
                                                                const List & logtheta,
                                                                const IntegerMatrix & V,
                                                                const NumericMatrix & alpha) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  int p = V.ncol();
  // 
  NumericMatrix clusterloglikelihoods(n, p);
  std::fill(clusterloglikelihoods.begin(), clusterloglikelihoods.end(), NA_REAL);
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      for (int field = 0; field < p; field++){
        clusterloglikelihoods(icluster, field) = ll_cluster_field(clsize[icluster],  clmembers(icluster,_),
                              theta[field], logtheta[field], V(_,field), alpha(icluster, field));
      }
    }
  }
  return clusterloglikelihoods;
}

// function to compute log likelihoods associated with each cluster and one particular field
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with NA
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p is "theta" in the paper
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// mVec which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_loglikelihood_all_clusters_one_field_cpp(int l ,
                                                               const List & clustering,
                                                               const NumericVector & theta_l,
                                                               const NumericVector & logtheta_l,
                                                               const IntegerMatrix & V,
                                                               const NumericMatrix & alpha) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  // 
  NumericVector cl_likelihood_field(n);
  std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), NA_REAL);
  
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      cl_likelihood_field(icluster) = ll_cluster_field(clsize[icluster],  clmembers(icluster,_),
                          theta_l, logtheta_l, V(_,l), alpha(icluster, l));
    }
  }
  return cl_likelihood_field;
}


