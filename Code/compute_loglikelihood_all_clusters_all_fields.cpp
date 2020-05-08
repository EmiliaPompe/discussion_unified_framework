#include <Rcpp.h>
using namespace Rcpp;

// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with NA
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericMatrix compute_loglikelihood_all_clusters_all_fields_cpp(const List & clustering,
                                             const List & theta,
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
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      NumericVector cl_likelihood_field(p);
      std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), 0.0);
      // compute likelihood recursively over members of cluster
      // first member (associated row in V)
      int j = clmembers(icluster,0);
      for (int l = 0; l < p; l++){
        NumericVector theta_field = theta[l];
        // recall p[cumdime[l] + V(j,l)] is theta_{l, v(j,l)} in the paper
        cl_likelihood_field(l) += log(theta_field(V(j,l)));
      }
      if (clsize[icluster] > 1){
        // if more members, implement recursion
        // loop over remaining members
        for (int imember = 1; imember < clsize[icluster]; imember++){
          // original index of that member (i.e. corresponding row in V)
          int q = clmembers(icluster,imember);
          // loop over fields 
          for (int l = 0; l < p; l++){
            NumericVector theta_field = theta[l];
            double logprod = 0.;
            // to implement first part of recursion
            // multiply by associated alpha' * V(q,l)
            cl_likelihood_field(l) += log(alpha(icluster,l) * theta_field(V(q,l)));
            // next, implement second part of recursion
            for (int othermember = 0; othermember < imember; othermember ++){
              int qprime = clmembers(icluster,othermember);
              int sameV = 0;
              if (V(q,l) == V(qprime,l)){
                sameV = 1;
              }
              logprod += log((1 - alpha(icluster,l)) * sameV + alpha(icluster,l) * theta_field(V(qprime,l)));
            }
            logprod += log((1 - alpha(icluster,l)) * theta_field(V(q,l)));
            // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
            double max_logs = std::max(logprod, cl_likelihood_field(l));
            cl_likelihood_field(l) = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field(l) - max_logs));
          }
        }
      }
      for (int l = 0; l < p; l++){
        clusterloglikelihoods(icluster,l) = cl_likelihood_field(l);
      }
    }
  }
  return clusterloglikelihoods;
}

// function to update log likelihoods associated with each cluster and *one single field* 
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with -Inf
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_loglikelihood_onefield_cpp(int field, 
                                    const List & clustering,
                                    const NumericVector & theta_field,
                                    const IntegerMatrix & V,
                                    const IntegerVector & Mvec,
                                    const NumericMatrix & alpha) {
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  NumericVector clusterloglikelihoods_onefield(n);
  std::fill(clusterloglikelihoods_onefield.begin(), clusterloglikelihoods_onefield.end(), R_NegInf);
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      double cl_likelihood_field = 0.;
      // compute likelihood recursively over members of cluster
      // first member (associated row in V)
      int j = clmembers(icluster,0);
      // recall p[cumdime[l] + V(j,l)] is theta_{l, v(j,l)} in the paper
      cl_likelihood_field += log(theta_field[V(j,field)]);
      if (clsize[icluster] > 1){
        // if more members, implement recursion
        // loop over remaining members
        for (int imember = 1; imember < clsize[icluster]; imember++){
          // original index of that member (i.e. corresponding row in V)
          int q = clmembers(icluster,imember);
          double logprod = 0.;
          // to implement first part of recursion
          // multiply by associated alpha' * V(q,l)
          cl_likelihood_field += log(alpha(icluster,field) * theta_field[V(q,field)]);
          // next, implement second part of recursion
          for (int othermember = 0; othermember < imember; othermember ++){
            int qprime = clmembers(icluster,othermember);
            int sameV = 0;
            if (V(q,field) == V(qprime,field)){
              sameV = 1;
            }
            logprod += log((1 - alpha(icluster,field)) * sameV + alpha(icluster,field) * theta_field[V(qprime,field)]);
          }
          logprod += log((1 - alpha(icluster,field)) * theta_field[V(q,field)]);
          // next we need to define exp(logprod) + exp(cl_likelihood_field)
          double max_logs = std::max(logprod, cl_likelihood_field);
          cl_likelihood_field = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field - max_logs));
        }
      }
      clusterloglikelihoods_onefield(icluster) = cl_likelihood_field;
    }
  }
  return clusterloglikelihoods_onefield;
}

