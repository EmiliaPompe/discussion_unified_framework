#ifndef _INCL_CMPTLL_
#define _INCL_CMPTLL_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

double compute_loglikelihood_one_cluster_one_field_cpp(int l ,
                                                       int icluster, 
                                                       const List & clustering,
                                                       const NumericVector & theta_l,
                                                       const NumericVector & logtheta_l,
                                                       const IntegerMatrix & V,
                                                       const double a);
  
NumericMatrix compute_loglikelihood_all_clusters_all_fields_cpp(const List & clustering,
                                                                const List & theta,
                                                                const List & logtheta,
                                                                const IntegerMatrix & V,
                                                                const NumericMatrix & alpha);
  
NumericVector compute_loglikelihood_all_clusters_one_field_cpp(int l ,
                                                               const List & clustering,
                                                               const NumericVector & theta_l,
                                                               const NumericVector & logtheta_l,
                                                               const IntegerMatrix & V,
                                                               const NumericMatrix & alpha);
  
#endif

