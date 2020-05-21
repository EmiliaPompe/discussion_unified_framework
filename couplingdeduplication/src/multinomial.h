#ifndef _INCL_MULTINOM_
#define _INCL_MULTINOM_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// takes a vector of log weights, a uniform draw
// and outputs a draw from the Categorical distribution with probabilities equal to the weights
int multinomial_(const NumericVector & logw, double uniform);

// takes two vectors of log weights, two uniform draws
// and outputs a pair of integers from a maximal coupling of Categorical distributions
IntegerVector coupled_multinomial_(const NumericVector & logw1, const NumericVector & logw2, double uniform1, double uniform2);

#endif

