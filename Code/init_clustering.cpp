#include <Rcpp.h>
using namespace Rcpp;

// Function to obtain partition/clustering based on vector lambda
// lambda must be a vector of length n with entries taking value in {0,...,n}

//  clustering returns 
//  * clize: a list of n values representing cluster sizes
//  * ksize: a number of non-empty clusters
//  * clmembers: an nxn matrix, where each row corresponds to a cluster
//    and each row 'i' contains, from left to right, the indices which(lambda == i)
//    (with C convention, so indices start at zero)
//    ... the rest of the matrix 'clmembers' is filled with '-1'

// [[Rcpp::export]]
List init_clustering(const IntegerVector & lambda) {
  int n = lambda.size();
  // int H = V.ncol();
  // initialize pclusterM, probM and psampq;
  int ksize = 0;
  IntegerVector clsize (n, 0);
  IntegerMatrix clmembers(n,n);
  std::fill(clmembers.begin(), clmembers.end(), -1);
  // compute clsize and ksize 
  for (int i = 0; i < n ; i++){
    int lambdai = lambda[i];
    if (clsize[lambdai] == 0){ 
      ksize += 1;
    }
    clmembers(lambdai, clsize[lambdai]) = i;
    clsize[lambdai] += 1;
  }
  return List::create(Named("clsize") = clsize, Named("ksize") = ksize, Named("clmembers") = clmembers);
}

