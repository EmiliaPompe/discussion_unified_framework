#ifndef _INCL_CLUSTERMANIP_
#define _INCL_CLUSTERMANIP_
#include <RcppEigen.h>
#include "ll_cluster_field.h"
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
                                 const NumericMatrix & alpha);

void compute_proba_eta(NumericMatrix & uponetajoining_loglikelihood, 
                        NumericVector & logproba_eta, 
                        NumericVector & newblock_loglikelihood, NumericVector & theta_field, NumericVector & logtheta_field, 
                        int ieta, int p, int n, 
                        const IntegerVector & clsize,
                        const IntegerMatrix & clmembers,
                        const NumericMatrix & clusterloglikelihoods,
                        const List & theta, const List & logtheta,
                        const IntegerMatrix & V,
                        const NumericMatrix & alpha,
                        const int  N, const int ksize);

#endif

