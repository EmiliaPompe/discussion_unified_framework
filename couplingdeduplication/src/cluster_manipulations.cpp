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