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
      remove_label_from_partition(ieta, label1,
                                  ksize1, clmembers1, clsize1,
                                  clusterloglikelihoods1,
                                  theta1, logtheta1, V, alpha1);
      remove_label_from_partition(ieta, label2,
                                  ksize2, clmembers2, clsize2,
                                  clusterloglikelihoods2,
                                  theta2, logtheta2, V, alpha2);

      // compute probability for new eta given the other variables 
      compute_proba_eta(uponetajoining_loglikelihood1, 
                        logproba_eta1, 
                        newblock_loglikelihood1, theta_field1, logtheta_field1, 
                        ieta, p, n, 
                        clsize1,
                        clmembers1,
                        clusterloglikelihoods1,
                        theta1, logtheta1,
                        V,
                        alpha1,
                        N1, ksize1);
      
      // compute probability for new eta given the other variables 
      compute_proba_eta(uponetajoining_loglikelihood2, 
                        logproba_eta2, 
                        newblock_loglikelihood2, theta_field2, logtheta_field2, 
                        ieta, p, n, 
                        clsize2,
                        clmembers2,
                        clusterloglikelihoods2,
                        theta2, logtheta2,
                        V,
                        alpha2,
                        N2, ksize2);
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
