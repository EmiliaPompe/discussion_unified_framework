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

// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times H = ncol(V)
// if a cluster is empty, the corresponding row is filled with zeros
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p which is "theta" in the paper
// V which is the observations
// dimV which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_loglikelihood_clusters(const List & clustering,
                                           const NumericVector & p,
                                           const IntegerMatrix & V,
                                           const IntegerVector & dimV,
                                           const NumericMatrix & a) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  int H = V.ncol();
  // 
  int cumdime[H]; 
  cumdime[0] = 0;
  for (int l=1; l<H; l++){
    cumdime[l] = cumdime[l-1] + dimV[l-1];
  }  
  //  
  NumericMatrix clusterloglikelihoods(n, H);
  std::fill(clusterloglikelihoods.begin(), clusterloglikelihoods.end(), 0.0);
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      NumericVector cl_likelihood_field(H);
      std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), 0.0);
      // compute likelihood recursively over members of cluster
      // first member (associated row in V)
      int j = clmembers[icluster,0];
      for (int l = 0; l < H; l++){
        // recall p[cumdime[l] + V(j,l)] is theta_{l, v(j,l)} in the paper
        cl_likelihood_field(l) += log(p[cumdime[l] + V(j,l)]);
      }
      if (clsize[icluster] > 1){
        // if more members, implement recursion
        // loop over remaining members
        for (int imember = 1; imember < clsize[icluster]; imember++){
          // original index of that member (i.e. corresponding row in V)
          int q = clmembers[icluster,imember];
          double logprod = 0.;
          // loop over fields 
          for (int l = 0; l < H; l++){
            // to implement first part of recursion
            // multiply by associated alpha' * V(q,l)
            cl_likelihood_field(l) += log(a(q,l) * p[cumdime[l] + V(q,l)]);
            // next, implement second part of recursion
            for (int othermember = 0; othermember < imember; othermember ++){
              int qprime = clmembers[icluster,othermember];
              logprod += log((1 - a(q,l)) * (V(q,l) == V(qprime,l)) + a(q,l) * p[cumdime[l] + V(qprime,l)]);
            }
            logprod += log((1 - a(q,l)) * p[cumdime[l] + V(q,l)]);
            // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
            double max_logs = std::max(logprod, cl_likelihood_field(l));
            cl_likelihood_field(l) = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field(l) - max_logs));
          }
        }
      }
      for (int l = 0; l < H; l++){
        clusterloglikelihoods(icluster,l) = cl_likelihood_field(l);
      }
    }
  }
  return clusterloglikelihoods;
}

// psampq is the ratio of include j vs exclude j
//    for (int l = 0 ; l < H; l ++){
//      probM[ l ] = 1.0;
//      pclusterM[ l ] = 1.0;
//      cumprod[ l ] = 1.0;
//    }
// pclusterM is the probability of excluding j in q
// probM is the probability of including j in q
//    bool firstElement = true;
//    for (int i = 0; i < n ; i++){
//      if (lambda[i] == q && j != i){
//        if (firstElement){
//          firstElement = false;
//          for (int l = 0; l < H; l++){
//            pclusterM[l] = p[cumdime[l] + V(j,l)];
//          }
/*        }else{
 for(int k = 0; k < i; k++){
 // for all the inspected record 
 if (lambda[k] == q && k != j){
 for(int l = 0 ; l < H; l++ ){
 cumprod[l] *= (1 - a(q,l)) * (V(k,l)== V(k,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
 pclusterM[l] = pclusterM[l] * a(q,l) * p[cumdime[l] + V(i,l)]+ (1 - a(q,l)) * p[cumdime[l] + V(i,l)] * cumprod[l];	
 }
 }
 }
}
 for (int l = 0; l < H; l++){
 probM[l] *= (1 - a(q,l)) * (V(i,l) == V(j,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
 }
}
}
 for (int l = 0; l < H; l++){
 probM[l] = p[cumdime[l] + V(j,l)] * ( a(q,l) + (1 - a(q,l)) * probM[l] / pclusterM[l]);
 }
 psampq[q] = 1.0;
 for (int l = 0; l < H; l++){
 psampq[q] *= probM[l];
 }
 }
 }
 
 return clusterlikelihoods;
 }
 */

/*  
 for (int q = 0; q < n; q++){
 // first see if q is a new cluster
 if (clsize[q] == 0){
 for (int l = 0; l < H; l++){
 psampq[q] *= p[cumdime[l]+ V(j,l)]; 
 }
 psampq[q] = psampq[q] * (N - ksize) / (n - ksize);
 }
 else{
 // if q is an observed cluster 
 // psampq is the ratio of include j vs exclude j
 for (int l = 0 ; l < H; l ++){
 probM[ l ] = 1.0;
 pclusterM[ l ] = 1.0;
 cumprod[ l ] = 1.0;
 }
 // pclusterM is the probability of excluding j in q
 // probM is the probability of including j in q
 bool firstElement = true;
 for (int i = 0; i < n ; i++){
 if (lambda[i] == q && j != i){
 if (firstElement){
 firstElement = false;
 for (int l = 0; l < H; l++){
 pclusterM[l] = p[cumdime[l] + V(j,l)];
 }
 }else{
 for(int k = 0; k < i; k++){
 // for all the inspected record 
 if (lambda[k] == q && k != j){
 for(int l = 0 ; l < H; l++ ){
 cumprod[l] *= (1 - a(q,l)) * (V(k,l)== V(k,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
 pclusterM[l] = pclusterM[l] * a(q,l) * p[cumdime[l] + V(i,l)]+ (1 - a(q,l)) * p[cumdime[l] + V(i,l)] * cumprod[l];	
 }
 }
 }
 }
 for (int l = 0; l < H; l++){
 probM[l] *= (1 - a(q,l)) * (V(i,l) == V(j,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
 }
 }
 }
 for (int l = 0; l < H; l++){
 probM[l] = p[cumdime[l] + V(j,l)] * ( a(q,l) + (1 - a(q,l)) * probM[l] / pclusterM[l]);
 }
 psampq[q] = 1.0;
 for (int l = 0; l < H; l++){
 psampq[q] *= probM[l];
 }
 }
 }
 
 // normalize psampq
 double sumpsampq = 0; 
 for (int q = 0; q < n ; q++){
 sumpsampq += psampq[q];
 }
 for (int q = 0; q < n ; q++){
 psampq[q] = psampq[q] / sumpsampq;
 }
 return clustering;
 }
 
 */
