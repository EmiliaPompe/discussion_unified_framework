#include <RcppEigen.h>
#include "ll_cluster_field.h"
using namespace Rcpp;
using namespace std;

// function that actually does the recursion
double ll_cluster_field(int clustersize, 
                        const IntegerVector & clmembers,
                        const NumericVector & theta_l, 
                        const NumericVector & logtheta_l, 
                        const IntegerVector & Vfield, const double a){
  double log_likelihood = NA_REAL ; 
  double logprod,max_logs;
  int q,j,qprime,sameV;
  double loga = log(a);
  double log1minusa = log(1-a);
  if (clustersize == 0){
    // do nothing and return NA
  } else {
    // compute likelihood recursively over members of cluster
    log_likelihood = 0.0;
    // first member (associated row in V)
    j = clmembers(0);
    // recall theta_l[ V(j,field)] is theta_{l, v(j,l)} in the paper
    log_likelihood += logtheta_l[Vfield(j)];
    if (log_likelihood > 1){
      // if more members, implement recursion
      // loop over remaining members
      for (int imember = 1; imember < clustersize; imember++){
        // original index of that member (i.e. corresponding row in V)
        q = clmembers(imember);
        logprod = 0.;
        // to implement first part of recursion
        // multiply by associated alpha' * theta_l[V(q,l)]
        log_likelihood += loga + logtheta_l[Vfield(q)];
        // next, implement second part of recursion
        for (int othermember = 0; othermember < imember; othermember ++){
          qprime = clmembers(othermember);
          sameV = 0;
          if (Vfield(q) == Vfield(qprime)){
            sameV = 1;
          }
          logprod += log((1 - a) * sameV + a * theta_l[Vfield(qprime)]);
        }
        logprod += log1minusa + logtheta_l[Vfield(q)];
        max_logs = std::max(logprod, log_likelihood);
        log_likelihood = max_logs + log(exp(logprod - max_logs) + exp(log_likelihood - max_logs));
      }
    }
  }
  return log_likelihood;  
}
