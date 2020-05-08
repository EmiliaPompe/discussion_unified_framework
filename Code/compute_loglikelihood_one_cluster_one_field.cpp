#include <Rcpp.h>
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
                                                   const IntegerMatrix & V,
                                                   const double a) {
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  // 
  double log_likelihood = NA_REAL ; 
  // if empty cluster, do nothing
  if (clsize[icluster] == 0){
    // do nothing and return 
    return log_likelihood;
  } else {
    // compute likelihood recursively over members of cluster
    log_likelihood = 0.0;
    // first member (associated row in V)
    int j = clmembers(icluster,0);
    // recall theta_l[ V(j,l)] is theta_{l, v(j,l)} in the paper
    log_likelihood += log(theta_l[ V(j,l)]);
    if (clsize[icluster] > 1){
      // if more members, implement recursion
      // loop over remaining members
      for (int imember = 1; imember < clsize[icluster]; imember++){
        // original index of that member (i.e. corresponding row in V)
        int q = clmembers(icluster,imember);
        // loop over fields 
        double logprod = 0.;
        // to implement first part of recursion
        // multiply by associated alpha' * theta_l[V(q,l)]
          log_likelihood += log(a * theta_l[ V(q,l)]);
          // next, implement second part of recursion
          for (int othermember = 0; othermember < imember; othermember ++){
            int qprime = clmembers(icluster,othermember);
            int sameV = 0;
            if (V(q,l) == V(qprime,l)){
              sameV = 1;
            }
            logprod += log((1 - a) * sameV + a * theta_l[ V(qprime,l)]);
          }
          logprod += log((1 - a) * theta_l[ V(q,l)]);
          // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
          // Rcout << "adds second term to l-th term" << logprod <<"\n";
          double max_logs = std::max(logprod, log_likelihood);
          log_likelihood = max_logs + log(exp(logprod - max_logs) + exp(log_likelihood - max_logs));
      }
    }
    return log_likelihood;  
  }
}

