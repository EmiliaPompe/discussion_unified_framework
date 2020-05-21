#include <RcppEigen.h>
#include "multinomial.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int multinomial_(const NumericVector & logw, double uniform){
  int n = logw.size();
  // exponentiate log weights, in stable way 
  double maxlogw = Rcpp::max(logw);
  NumericVector cumsum_w = exp(logw - maxlogw);
  int draw = 0;
  // cumsum_w = wrap(cumsum(logw));
  for (int i = 1; i < n; i++){
    cumsum_w(i) = cumsum_w(i) + cumsum_w(i-1);
  }
  uniform = uniform * cumsum_w(n-1);
  int running_index = 0;
  double sumw = cumsum_w(0);
  if (uniform <= sumw){
    draw = running_index;
  } else {
    while (uniform > sumw){
      sumw = cumsum_w(running_index);
      running_index ++;
    }
    draw = running_index-1;
  }
  return draw;
}

// [[Rcpp::export]]
IntegerVector coupled_multinomial_(const NumericVector & logw1, const NumericVector & logw2, double uniform1, double uniform2){
  int n = logw1.size();
  IntegerVector draws(2);
  // exponentiate log weights, in stable way 
  double maxlogw1 = Rcpp::max(logw1);
  NumericVector w1 = exp(logw1 - maxlogw1);
  w1 = w1 / sum(w1);
  double maxlogw2 = Rcpp::max(logw2);
  NumericVector w2 = exp(logw2 - maxlogw2);
  w2 = w2 / sum(w2);
  // compute common part between two distributions
  NumericVector commonpart = Rcpp::pmin(w1, w2);
  double sumcommonpart = sum(commonpart);
  // sumcommonpart = 1 - TV between two distributions
  if (uniform1 < sumcommonpart){
    // draw from common part
    NumericVector cumsumw = cumsum(commonpart);
    double u = uniform2 * cumsumw(n-1);
    int running_index = 0;
    double sumw = cumsumw(running_index);
    if (u <= sumw){
      draws[0] = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw(running_index);
        running_index ++;
      }
      draws[0] = running_index-1;
    }
    // set second draw equal to first draw
    draws[1] = draws[0];
  } else {
    // otherwise draw from residuals, w1 - commonpart, and w2 - commonpart
    NumericVector cumsumw1 = cumsum(w1 - commonpart);
    double u = uniform2 * cumsumw1(n-1);
    int running_index = 0;
    double sumw = cumsumw1(0);
    if (u <= sumw){
      draws[0] = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw1(running_index);
        running_index ++;
      }
      draws[0] = running_index-1;
    }
    NumericVector cumsumw2 = cumsum(w2 - commonpart);
    u = uniform2 * cumsumw2(n-1);
    running_index = 0;
    sumw = cumsumw2(0);
    if (u <= sumw){
      draws[1] = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw2(running_index);
        running_index ++;
      }
      draws[1] = running_index-1;
    }
  }
  return draws;
}
