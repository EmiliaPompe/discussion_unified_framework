#include <Rcpp.h>
using namespace Rcpp;

// Function to obtain partition/clustering based on vector eta
// eta must be a vector of length n with entries taking value in {0,...,n-1}

//  clustering returns 
//  * clize: a list of n values representing cluster sizes
//  * ksize: a number of non-empty clusters
//  * clmembers: an nxn matrix, where each row corresponds to a cluster
//    and each row 'i' contains, from left to right, the indices which(eta == i)
//    (with C convention, so indices start at zero)
//    ... the rest of the matrix 'clmembers' is filled with '-1'

// [[Rcpp::export]]
List init_clustering(const IntegerVector & eta) {
  int n = eta.size();
  int ksize = 0;
  IntegerVector clsize (n, 0);
  IntegerMatrix clmembers(n,n);
  std::fill(clmembers.begin(), clmembers.end(), -1);
  // compute clsize and ksize 
  for (int indiv = 0; indiv < n ; indiv++){
    int eta_indiv = eta[indiv];
    if (clsize[eta_indiv] == 0){ 
      ksize += 1;
    }
    clmembers(eta_indiv, clsize[eta_indiv]) = indiv;
    clsize[eta_indiv] += 1;
  }
  return List::create(Named("clsize") = clsize, Named("ksize") = ksize, Named("clmembers") = clmembers);
}



