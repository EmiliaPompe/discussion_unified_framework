#ifndef _INCL_LLCLUSTERFIELD_
#define _INCL_LLCLUSTERFIELD_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// compute log-likelihood associated with one cluster and one field
double ll_cluster_field(int clustersize, 
                        const IntegerVector & clmembers,
                        const NumericVector & theta_l, 
                        const NumericVector & logtheta_l, 
                        const IntegerVector & Vfield, const double a);
#endif

