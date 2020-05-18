#include <Rcpp.h>
#include "ll_cluster_field.h"
#include "compute_loglikelihood.h"
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
List update_eta_cpp(IntegerVector & eta,
                    List & clustering,
                    const NumericMatrix previous_clusterloglikelihoods,
                    const List & theta,
                    const List & logtheta,
                    const IntegerMatrix & V,
                    const NumericMatrix & alpha,
                    int N) {
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  IntegerVector ksize_vec = clustering["ksize"];
  int ksize = ksize_vec(0);
  //
  int n = V.nrow();
  int p = V.ncol();
  NumericMatrix clusterloglikelihoods = clone(previous_clusterloglikelihoods);
  // 
  NumericVector newblock_loglikelihood = no_init(p);
  NumericVector theta_field;
  NumericVector logtheta_field;
  NumericMatrix uponetajoining_loglikelihood = no_init(n,p);
  NumericVector logproba_eta = no_init(n);
  NumericVector uniform;
  
  // loop over components of eta
  for (int ieta = 0; ieta < n; ieta++){
    int label = eta[ieta];
    // first, remove eta[ieta] from current partition
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
      // Rcout << "putting NA in row" << label << std::endl;
      // decrement cluster counter
      ksize --;
      // remove log-likelihoods associated with that cluster
      std::fill(clusterloglikelihoods.row(label).begin(), clusterloglikelihoods.row(label).end(), NA_REAL);
      // for (int field = 0; field < p; field ++){
      //   clusterloglikelihoods(label,field) = NA_REAL;
      // }
    } else {
      // Rcout << "recomputing row " << label << std::endl;
      // recompute likelihood of cluster from scratch
      for (int field = 0; field < p; field ++){
        clusterloglikelihoods(label,field) = compute_loglikelihood_one_cluster_one_field_cpp(field,
                              label, clustering, theta[field], logtheta[field], V, alpha(label,field));
      }        
    }
    // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
    // next we can compute the likelihood associated with a new block with that label in it
    for (int l = 0; l < p; l++){
      logtheta_field = logtheta[l];
      newblock_loglikelihood(l) = logtheta_field(V(ieta,l));
    }
    // vector to aggregate likelihood and prior probabilities for new label
    std::fill(logproba_eta.begin(), logproba_eta.end(), 0.);
    // next we loop over clusters
    for (int icluster = 0; icluster < n; icluster ++){
      // first compute likelihoods of joining
      if (clsize[icluster]==0){
        // this would be a new cluster
        uponetajoining_loglikelihood.row(icluster) = newblock_loglikelihood;
      } else {
        // this would be an existing cluster
        for (int l = 0; l < p; l++){
          theta_field = theta[l];
          logtheta_field = logtheta[l];
          // first part of recursion
          uponetajoining_loglikelihood(icluster, l) = log(alpha(icluster,l)) + logtheta_field(V(ieta,l)) + 
            clusterloglikelihoods(icluster,l);
          double logprod = 0.;
          // next, implement second part of recursion
          for (int clustermember = 0; clustermember < clsize[icluster]; clustermember ++){
            int qprime = clmembers(icluster, clustermember);
            logprod += log((1 - alpha(icluster,l)) * (V(ieta,l) == V(qprime,l)) + alpha(icluster,l) * theta_field(V(qprime,l)));
          }
          logprod += log(1 - alpha(icluster,l)) + logtheta_field(V(ieta,l));
          double max_logs = std::max(logprod, uponetajoining_loglikelihood(icluster, l));
          uponetajoining_loglikelihood(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood(icluster, l) - max_logs));
        }
      }
      // then aggregate priors and conditional likelihood
      if (clsize[icluster] == 0){
        // i.e. P(eta[ieta] = q) where q is the label of a new block
        logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster));
        // logproba_eta[icluster] += (log(N-ksize)); // CAREFUL CHECK - log(n-ksize));
        logproba_eta[icluster] += (log(N-ksize) - log(n-ksize)); // CAREFUL CHECK!!
        // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster 
      } else {
        logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster) - clusterloglikelihoods.row(icluster));
        logproba_eta[icluster] += 0.;
      }
    }
    // sample from categorical/multinomial given logproba_eta
    double maxlogproba_eta = Rcpp::max(logproba_eta);
    logproba_eta = exp(logproba_eta - maxlogproba_eta);
    int draw = 0;
    NumericVector cumsumw = cumsum(logproba_eta);
    GetRNGstate();
    uniform = runif(1);
    PutRNGstate();
    double u = uniform(0) * cumsumw(n-1);
    int running_index = 0;
    double sumw = cumsumw(running_index);
    if (u <= sumw){
      draw = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw(running_index);
        running_index ++;
      }
      draw = running_index-1;
    }
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
  clustering["clsize"] = clsize;
  clustering["ksize"] = ksize;
  clustering["clmembers"] = clmembers;
  return List::create(Named("eta") = eta, 
                      Named("clusterloglikelihoods") = clusterloglikelihoods, 
                      Named("clustering") = clustering);
}

