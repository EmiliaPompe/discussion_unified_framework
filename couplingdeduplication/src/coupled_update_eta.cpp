#include <Rcpp.h>
#include "ll_cluster_field.h"
#include "compute_loglikelihood.h"
#include "cluster_manipulations.h"
#include "multinomial.h"
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List coupled_update_eta(
    IntegerVector & eta1,
    IntegerVector & eta2,
    List & clustering1,
    List & clustering2,
    const NumericMatrix & previous_clusterloglikelihoods1,
    const NumericMatrix & previous_clusterloglikelihoods2,
    const List & theta1,
    const List & theta2,
    const List & logtheta1,
    const List & logtheta2,
    const IntegerMatrix & V,
    const NumericMatrix & alpha1,
    const NumericMatrix & alpha2,
    int N1, int N2, double updateprobability) {
  // Rcout << "start" << endl;
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
  //
  NumericVector newblock_loglikelihood1 = no_init(p);
  NumericVector newblock_loglikelihood2 = no_init(p);
  NumericVector theta_field1, theta_field2;
  NumericVector logtheta_field1, logtheta_field2;
  NumericMatrix uponetajoining_loglikelihood1 = no_init(n,p);
  NumericMatrix uponetajoining_loglikelihood2 = no_init(n,p);
  NumericVector logproba_eta1 = no_init(n);
  NumericVector logproba_eta2 = no_init(n);
  NumericVector firstuniform, seconduniform;
  // to decide whether to update or not
  GetRNGstate();
  NumericVector unifs_ = runif(n);
  PutRNGstate();
  // loop over components of eta
  for (int ieta = 0; ieta < n; ieta++){
    if (unifs_[ieta] < updateprobability){
      int label1 = eta1[ieta];
      int label2 = eta2[ieta];
      // Rcout << "removing labels..." << endl;
      remove_label_from_partition(ieta, label1,
                                  ksize1, clmembers1, clsize1,
                                  clusterloglikelihoods1,
                                  theta1, logtheta1, V, alpha1);
      // Rcout << "... partially done..." << endl;
      remove_label_from_partition(ieta, label2,
                                  ksize2, clmembers2, clsize2,
                                  clusterloglikelihoods2,
                                  theta2, logtheta2, V, alpha2);
      // Rcout << "... done" << endl;
      // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
      // next we can compute the likelihood associated with a new block with that label in it
      for (int field = 0; field < p; field++){
        logtheta_field1 = logtheta1[field];
        newblock_loglikelihood1(field) = logtheta_field1(V(ieta,field));
        logtheta_field2 = logtheta2[field];
        newblock_loglikelihood2(field) = logtheta_field2(V(ieta,field));
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
          for (int field = 0; field < p; field++){
            theta_field1 = theta1[field];
            logtheta_field1 = logtheta1[field];
            // first part of recursion
            uponetajoining_loglikelihood1(icluster, field) = log(alpha1(icluster,field)) + logtheta_field1(V(ieta,field)) +
              clusterloglikelihoods1(icluster,field);
            double logprod = 0.;
            // next, implement second part of recursion
            for (int clustermember = 0; clustermember < clsize1[icluster]; clustermember ++){
              int qprime = clmembers1(icluster, clustermember);
              logprod += log((1 - alpha1(icluster,field)) * (V(ieta,field) == V(qprime,field)) + alpha1(icluster,field) * theta_field1(V(qprime,field)));
            }
            logprod += log(1 - alpha1(icluster,field)) + logtheta_field1(V(ieta,field));
            double max_logs = std::max(logprod, uponetajoining_loglikelihood1(icluster, field));
            uponetajoining_loglikelihood1(icluster, field) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood1(icluster, field) - max_logs));
          }
        }
        // first compute likelihoods of joining
        if (clsize2[icluster]==0){
          // this would be a new cluster
          uponetajoining_loglikelihood2.row(icluster) = newblock_loglikelihood2;
        } else {
          // this would be an existing cluster
          for (int field = 0; field < p; field++){
            theta_field2 = theta2[field];
            logtheta_field2 = logtheta2[field];
            // first part of recursion
            uponetajoining_loglikelihood2(icluster, field) = log(alpha2(icluster,field)) + logtheta_field2(V(ieta,field)) +
              clusterloglikelihoods2(icluster,field);
            double logprod = 0.;
            // next, implement second part of recursion
            for (int clustermember = 0; clustermember < clsize2[icluster]; clustermember ++){
              int qprime = clmembers2(icluster, clustermember);
              logprod += log((1 - alpha2(icluster,field)) * (V(ieta,field) == V(qprime,field)) + alpha2(icluster,field) * theta_field2(V(qprime,field)));
            }
            logprod += log(1 - alpha2(icluster,field)) + logtheta_field2(V(ieta,field));
            double max_logs = std::max(logprod, uponetajoining_loglikelihood2(icluster, field));
            uponetajoining_loglikelihood2(icluster, field) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood2(icluster, field) - max_logs));
          }
        }
        // then aggregate priors and conditional likelihood
        if (clsize1[icluster] == 0){
          logproba_eta1[icluster] = sum(uponetajoining_loglikelihood1.row(icluster));
          logproba_eta1[icluster] += (log(N1-ksize1) - log(n-ksize1)); 
        } else {
          logproba_eta1[icluster] = sum(uponetajoining_loglikelihood1.row(icluster) - clusterloglikelihoods1.row(icluster));
          logproba_eta1[icluster] += 0.;
        }
        if (clsize2[icluster] == 0){
          logproba_eta2[icluster] = sum(uponetajoining_loglikelihood2.row(icluster));
          logproba_eta2[icluster] += (log(N2-ksize2) - log(n-ksize2)); 
        } else {
          logproba_eta2[icluster] = sum(uponetajoining_loglikelihood2.row(icluster) - clusterloglikelihoods2.row(icluster));
          logproba_eta2[icluster] += 0.;
        }
      }
      
      // sample from maximal coupling given logproba_eta1, logproba_eta2
      GetRNGstate();
      NumericVector uniforms = runif(2); 
      PutRNGstate();
      IntegerVector draws = coupled_multinomial_(logproba_eta1, logproba_eta2, uniforms(0), uniforms(1));
      // Rcout << draws << endl;
      // thus we update eta, and the clustering, and associated loglikelihoods appropriately
      eta1[ieta] = draws[0];
      eta2[ieta] = draws[1];
      // update cluster log likelihoods
      clusterloglikelihoods1.row(draws[0]) = uponetajoining_loglikelihood1.row(draws[0]);
      clusterloglikelihoods2.row(draws[1]) = uponetajoining_loglikelihood2.row(draws[1]);
      // if cluster is new, we need to draw a vector alpha
      // and we need to update the partition
      // change number of cluster if need be
      if (clsize1[draws[0]] == 0){
        ksize1 ++;
      }
      if (clsize2[draws[1]] == 0){
        ksize2 ++;
      }
      // adds index in matrix storing members of each cluster
      clmembers1(draws[0], clsize1[draws[0] ]) = ieta;
      clmembers2(draws[1], clsize2[draws[1] ]) = ieta;
      // increment cluster size
      clsize1[draws[0] ] += 1;
      clsize2[draws[1] ] += 1;
    }
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
