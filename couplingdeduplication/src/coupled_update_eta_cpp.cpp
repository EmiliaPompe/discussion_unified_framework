#include <Rcpp.h>
#include "ll_cluster_field.h"
#include "compute_loglikelihood.h"
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List coupled_update_eta_cpp(
  IntegerVector & eta1,
  IntegerVector & eta2,
  List & clustering1,
  List & clustering2,
  const NumericMatrix & previous_clusterloglikelihoods1,
  const NumericMatrix & previous_clusterloglikelihoods2,
  NumericMatrix & uponetajoining_loglikelihood1,
  NumericMatrix & uponetajoining_loglikelihood2,
  const List & theta1,
  const List & theta2,
  const List & logtheta1,
  const List & logtheta2,
  const IntegerMatrix & V,
  const NumericMatrix & alpha1,
  const NumericMatrix & alpha2,
  int N1, int N2) {
  IntegerVector clsize1    = clustering1["clsize"];
  IntegerMatrix clmembers1 = clustering1["clmembers"];
  IntegerVector ksize_vec1 = clustering1["ksize"];
  IntegerVector clsize2    = clustering2["clsize"];
  IntegerMatrix clmembers2 = clustering2["clmembers"];
  IntegerVector ksize_vec2 = clustering2["ksize"];
  int ksize1 = ksize_vec1(0);
  int ksize2 = ksize_vec2(0);
  int n = V.nrow();
  int p = V.ncol();
  //
    NumericMatrix clusterloglikelihoods1 = clone(wrap(previous_clusterloglikelihoods1));
    NumericMatrix clusterloglikelihoods2 = clone(wrap(previous_clusterloglikelihoods2));
    NumericVector newblock_loglikelihood1 = no_init(p);
    NumericVector newblock_loglikelihood2 = no_init(p);
    NumericVector theta_field1, theta_field2;
    NumericVector logtheta_field1, logtheta_field2;
    // NumericMatrix uponetajoining_loglikelihood = no_init(n,p);
    NumericVector logproba_eta1 = no_init(n);
    NumericVector logproba_eta2 = no_init(n);
    NumericVector firstuniform, seconduniform;
    
    // loop over components of eta
    for (int ieta = 0; ieta < n; ieta++){
      int label1 = eta1[ieta];
      int label2 = eta2[ieta];
      // first, remove eta[ieta] from current partition
      // i.e. translate -1's towards the left
    bool shift = false;
    for (int icol = 0; icol < clsize1[label1]; icol++){
      if (clmembers1(label1,icol) == ieta){
        shift = true;
      }
      if (shift){
        // move entry of next column to current column
        if (icol + 1 < n){
          clmembers1(label1,icol) = clmembers1(label1,icol+1);
        } else {
          clmembers1(label1,icol) = -1;
        }
      }
    }
    for (int icol = 0; icol < clsize2[label2]; icol++){
      if (clmembers2(label2,icol) == ieta){
        shift = true;
      }
      if (shift){
        // move entry of next column to current column
        if (icol + 1 < n){
          clmembers2(label2,icol) = clmembers2(label2,icol+1);
        } else {
          clmembers2(label2,icol) = -1;
        }
      }
    }

    clsize1[label1] --;
    clsize2[label2] --;
    // if cluster is now empty
    if (clsize1[label1] == 0){
      // decrement cluster counter
      ksize1 --;
      // remove log-likelihoods associated with that cluster
      std::fill(clusterloglikelihoods1.row(label1).begin(), clusterloglikelihoods1.row(label1).end(), NA_REAL);
      // remove corresponding alpha
    } else {
      // recompute likelihood of cluster from scratch
      for (int field = 0; field < p; field ++){
        clusterloglikelihoods1(label1,field) = compute_loglikelihood_one_cluster_one_field_cpp(field,
                              label1, clustering1, theta1[field], logtheta1[field], V, alpha1(label1,field));
      }
    }
    // if cluster is now empty
    if (clsize2[label2] == 0){
      // decrement cluster counter
      ksize2 --;
      // remove log-likelihoods associated with that cluster
      std::fill(clusterloglikelihoods2.row(label2).begin(), clusterloglikelihoods2.row(label2).end(), NA_REAL);
      // remove corresponding alpha
    } else {
      // recompute likelihood of cluster from scratch
      for (int field = 0; field < p; field ++){
        clusterloglikelihoods2(label2,field) = compute_loglikelihood_one_cluster_one_field_cpp(field,
                              label2, clustering2, theta2[field], logtheta2[field], V, alpha2(label2,field));
      }
    }
    // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
    // next we can compute the likelihood associated with a new block with that label in it
    for (int l = 0; l < p; l++){
      logtheta_field1 = logtheta1[l];
      newblock_loglikelihood1(l) = logtheta_field1(V(ieta,l));
      logtheta_field2 = logtheta2[l];
      newblock_loglikelihood2(l) = logtheta_field2(V(ieta,l));
    }
    // vector to aggregate likelihood and prior probabilities for new label
    std::fill(logproba_eta1.begin(), logproba_eta1.end(), 0.);
    std::fill(logproba_eta2.begin(), logproba_eta2.end(), 0.);
    // next we loop over clusters
    for (int icluster = 0; icluster < n; icluster ++){
      // first compute likelihoods of joining
      if (clsize1[icluster]==0){
        // this would be a new cluster
        uponetajoining_loglikelihood1.row(icluster) = newblock_loglikelihood1;
      } else {
        // this would be an existing cluster
        for (int l = 0; l < p; l++){
          theta_field1 = theta1[l];
          logtheta_field1 = logtheta1[l];
          // first part of recursion
          uponetajoining_loglikelihood1(icluster, l) = log(alpha1(icluster,l)) + logtheta_field1(V(ieta,l)) +
            clusterloglikelihoods1(icluster,l);
          double logprod = 0.;
          // next, implement second part of recursion
          for (int clustermember = 0; clustermember < clsize1[icluster]; clustermember ++){
            int qprime = clmembers1(icluster, clustermember);
            logprod += log((1 - alpha1(icluster,l)) * (V(ieta,l) == V(qprime,l)) + alpha1(icluster,l) * theta_field1(V(qprime,l)));
          }
          logprod += log(1 - alpha1(icluster,l)) + logtheta_field1(V(ieta,l));
          double max_logs = std::max(logprod, uponetajoining_loglikelihood1(icluster, l));
          uponetajoining_loglikelihood1(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood1(icluster, l) - max_logs));
        }
      }
      // first compute likelihoods of joining
      if (clsize2[icluster]==0){
        // this would be a new cluster
        uponetajoining_loglikelihood2.row(icluster) = newblock_loglikelihood2;
      } else {
        // this would be an existing cluster
        for (int l = 0; l < p; l++){
          theta_field2 = theta2[l];
          logtheta_field2 = logtheta2[l];
          // first part of recursion
          uponetajoining_loglikelihood2(icluster, l) = log(alpha2(icluster,l)) + logtheta_field2(V(ieta,l)) +
            clusterloglikelihoods2(icluster,l);
          double logprod = 0.;
          // next, implement second part of recursion
          for (int clustermember = 0; clustermember < clsize2[icluster]; clustermember ++){
            int qprime = clmembers2(icluster, clustermember);
            logprod += log((1 - alpha2(icluster,l)) * (V(ieta,l) == V(qprime,l)) + alpha2(icluster,l) * theta_field2(V(qprime,l)));
          }
          logprod += log(1 - alpha2(icluster,l)) + logtheta_field2(V(ieta,l));
          double max_logs = std::max(logprod, uponetajoining_loglikelihood2(icluster, l));
          uponetajoining_loglikelihood2(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood2(icluster, l) - max_logs));
        }
      }
      // then aggregate priors and conditional likelihood
      if (clsize1[icluster] == 0){
        // i.e. P(eta[ieta] = q) where q is the label of a new block
        logproba_eta1[icluster] = sum(uponetajoining_loglikelihood1.row(icluster));
        // logproba_eta[icluster] += (log(N-ksize)); // CAREFUL CHECK - log(n-ksize));
        logproba_eta1[icluster] += (log(N1-ksize1) - log(n-ksize1)); // CAREFUL CHECK!!
        // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster
      } else {
        logproba_eta1[icluster] = sum(uponetajoining_loglikelihood1.row(icluster) - clusterloglikelihoods1.row(icluster));
        logproba_eta1[icluster] += 0.;
      }
      // then aggregate priors and conditional likelihood
      if (clsize2[icluster] == 0){
        // i.e. P(eta[ieta] = q) where q is the label of a new block
        logproba_eta2[icluster] = sum(uponetajoining_loglikelihood2.row(icluster));
        // logproba_eta[icluster] += (log(N-ksize)); // CAREFUL CHECK - log(n-ksize));
        logproba_eta2[icluster] += (log(N2-ksize2) - log(n-ksize2)); // CAREFUL CHECK!!
        // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster
      } else {
        logproba_eta2[icluster] = sum(uponetajoining_loglikelihood2.row(icluster) - clusterloglikelihoods2.row(icluster));
        logproba_eta2[icluster] += 0.;
      }
    }

    // sample from maximal coupling given logproba_eta1, logproba_eta2
    double maxlogproba_eta1 = Rcpp::max(logproba_eta1);
    double maxlogproba_eta2 = Rcpp::max(logproba_eta2);
    logproba_eta1 = exp(logproba_eta1 - maxlogproba_eta1);
    logproba_eta2 = exp(logproba_eta2 - maxlogproba_eta2);
    // compute common part
    NumericVector commonpart = Rcpp::pmin(logproba_eta1, logproba_eta2);
    double sumcommonpart = sum(commonpart);
    GetRNGstate();
    firstuniform = runif(1);
    seconduniform = runif(1);
    PutRNGstate();
    int draw1 = 0;
    int draw2 = 0;
    if (firstuniform(0) < sumcommonpart){
      // draw from common part
      NumericVector cumsumw = cumsum(commonpart);
      double u = seconduniform(0) * cumsumw(n-1);
      int running_index = 0;
      double sumw = cumsumw(running_index);
      if (u <= sumw){
        draw1 = running_index;
      } else {
        while (u > sumw){
          sumw = cumsumw(running_index);
          running_index ++;
        }
        draw1 = running_index-1;
      }
      draw2 = draw1;
    } else {
      // draw from residuals
      NumericVector cumsumw1 = cumsum(logproba_eta1 - commonpart);
      double u = seconduniform(0) * cumsumw1(n-1);
      int running_index = 0;
      double sumw = cumsumw1(running_index);
      if (u <= sumw){
        draw1 = running_index;
      } else {
        while (u > sumw){
          sumw = cumsumw1(running_index);
          running_index ++;
        }
        draw1 = running_index-1;
      }
      NumericVector cumsumw2 = cumsum(logproba_eta2 - commonpart);
      u = seconduniform(0) * cumsumw2(n-1);
      running_index = 0;
      sumw = cumsumw2(running_index);
      if (u <= sumw){
        draw2 = running_index;
      } else {
        while (u > sumw){
          sumw = cumsumw2(running_index);
          running_index ++;
        }
        draw2 = running_index-1;
      }
    }
    // thus we update eta, and the clustering, and associated loglikelihoods appropriately
    eta1[ieta] = draw1;
    eta2[ieta] = draw2;
    // update cluster log likelihoods
    clusterloglikelihoods1.row(draw1) = uponetajoining_loglikelihood1.row(draw1);
    clusterloglikelihoods2.row(draw2) = uponetajoining_loglikelihood2.row(draw2);
    // if cluster is new, we need to draw a vector alpha
    // and we need to update the partition
    // change number of cluster if need be
    if (clsize1[draw1] == 0){
      ksize1 ++;
    }
    if (clsize2[draw2] == 0){
      ksize2 ++;
    }
    // adds index in matrix storing members of each cluster
    clmembers1(draw1, clsize1[draw1]) = ieta;
    clmembers2(draw2, clsize2[draw2]) = ieta;
    // increment cluster size
    clsize1[draw1] += 1;
    clsize2[draw2] += 1;
  }
  clustering1["clsize"]    = clsize1;
  clustering1["ksize"]     = ksize1;
  clustering1["clmembers"] = clmembers1;
  clustering2["clsize"]    = clsize2;
  clustering2["ksize"]     = ksize2;
  clustering2["clmembers"] = clmembers2;
  return List::create(Named("eta1") = eta1, 
                      Named("eta2") = eta2,
                      Named("clusterloglikelihoods1") = clusterloglikelihoods1,
                      Named("clusterloglikelihoods2") = clusterloglikelihoods2,
                      Named("clustering1") = clustering1,
                      Named("clustering2") = clustering2);
}
