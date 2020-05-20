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
#endif

