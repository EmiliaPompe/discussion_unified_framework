#include <Rcpp.h>
#include "ll_cluster_field.h"
#include "compute_loglikelihood.h"
#include "cluster_manipulations.h"
#include "multinomial.h"
using namespace Rcpp;
using namespace std;

// function to update the eta vector
// where eta[ieta] is a label referring to one of N clusters
// takes as arguments 
// * eta_previous, the vector of current etas
// * clustering_previous, the associated clustering
// * ll_previous, the log-likelihoods associated with each cluster
// * theta, list of frequencies for each field
// * logtheta, list of log-frequencies for each field
// * V, observations
// * alpha, corruption probabilities for each field and cluster
// * N
// [[Rcpp::export]]
List update_eta(IntegerVector & eta,
                List & clustering,
                const NumericMatrix previous_clusterloglikelihoods,
                const List & theta,
                const List & logtheta,
                const IntegerMatrix & V,
                const NumericMatrix & alpha,
                int N, double updateprobability) {
  IntegerVector clsize    = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  IntegerVector ksize_vec = clustering["ksize"];
  int ksize = ksize_vec(0);
  //
  int n = V.nrow();
  int p = V.ncol();
  NumericMatrix clusterloglikelihoods = clone(previous_clusterloglikelihoods);
  // 
  NumericVector newblock_loglikelihood = no_init(p);
  NumericMatrix uponetajoining_loglikelihood = no_init(n,p);
  NumericVector logproba_eta = no_init(n);
  NumericVector theta_field;
  NumericVector logtheta_field;
  NumericVector uniform;
  GetRNGstate();
  NumericVector unifs_ = runif(n);
  PutRNGstate();
  // loop over components of eta
  for (int ieta = 0; ieta < n; ieta++){
    if (unifs_[ieta] < updateprobability){
      int label = eta[ieta];
      // first, remove eta[ieta] from current partition
      remove_label_from_partition(ieta, label,
                                  ksize, clmembers, clsize,
                                  clusterloglikelihoods,
                                  theta, logtheta, V, alpha);
      // compute probability for new eta given the other variables 
      compute_proba_eta(uponetajoining_loglikelihood, 
                        logproba_eta, 
                        newblock_loglikelihood, theta_field, logtheta_field, 
                        ieta, p, n, 
                        clsize,
                        clmembers,
                        clusterloglikelihoods,
                        theta, logtheta,
                        V,
                        alpha,
                        N, ksize);
      // sample from categorical/multinomial given logproba_eta
      GetRNGstate();
      NumericVector u = runif(1);
      PutRNGstate();
      int draw = multinomial_(logproba_eta, u(0));
      // thus we update eta, and the clustering, and associated loglikelihoods appropriately
      eta[ieta] = draw;
      // update cluster log likelihoods
      // Rcout << "filling row " << draw << std::endl;
      clusterloglikelihoods.row(draw) = uponetajoining_loglikelihood.row(draw);
      // if cluster is new, we need to draw a vector alpha 
      // and we need to update the partition
      // change number of cluster if need be
      if (clsize[draw] == 0){ 
        ksize += 1;
      }
      // adds index in matrix storing members of each cluster
      clmembers(draw, clsize[draw]) = ieta;
      // increment cluster size
      clsize[draw] += 1;
    }
  }
  clustering["clsize"] = clsize;
  clustering["ksize"] = ksize;
  clustering["clmembers"] = clmembers;
  return List::create(Named("eta") = eta, 
                      Named("clusterloglikelihoods") = clusterloglikelihoods, 
                      Named("clustering") = clustering);
}

