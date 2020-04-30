#include <Rcpp.h>
using namespace Rcpp;

// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times H = ncol(V)
// if a cluster is empty, the corresponding row is filled with zeros
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p is "theta" in the paper
// V is the observations (careful, values need to start at zero and go to values in dimV-1)
// dimV which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_loglikelihood_clusters_field(int l ,
                                            const List & clustering,
                                             const NumericVector & thetaL,
                                             const IntegerMatrix & V,
                                             const IntegerVector & dimV,
                                             const NumericMatrix & alpha) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  // 
  NumericVector cl_likelihood_field(n);
  std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), 0.0);
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      // compute likelihood recursively over members of cluster
      // first member (associated row in V)
      int j = clmembers(icluster,0);
      // recall thetaL[ V(j,l)] is theta_{l, v(j,l)} in the paper
      cl_likelihood_field[icluster] += log(thetaL[ V(j,l)]);
      if (clsize[icluster] > 1){
        // if more members, implement recursion
        // loop over remaining members
        for (int imember = 1; imember < clsize[icluster]; imember++){
          // original index of that member (i.e. corresponding row in V)
          int q = clmembers(icluster,imember);
          // loop over fields 
            double logprod = 0.;
            // to implement first part of recursion
            // multiply by associated alpha' * V(q,l)
            cl_likelihood_field[icluster] += log(alpha(icluster,l) * thetaL[ V(q,l)]);
            // next, implement second part of recursion
            for (int othermember = 0; othermember < imember; othermember ++){
              int qprime = clmembers(icluster,othermember);
              int sameV = 0;
              if (V(q,l) == V(qprime,l)){
                sameV = 1;
              }
              logprod += log((1 - alpha(icluster,l)) * sameV + alpha(icluster,l) * thetaL[ V(qprime,l)]);
            }
            logprod += log((1 - alpha(icluster,l)) * thetaL[ V(q,l)]);
            // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
            // Rcout << "adds second term to l-th term" << logprod <<"\n";
            double max_logs = std::max(logprod, cl_likelihood_field[icluster]);
            cl_likelihood_field[icluster] = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field[icluster] - max_logs));
        }
      }
    }
  }
  return cl_likelihood_field;
}

