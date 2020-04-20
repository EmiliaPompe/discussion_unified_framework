#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix computeNweights(double g,
	int N1,
	int N2,
	int stepsize,
	int ksize1,
	int ksize2,
	int n) {
	// Rprintf("N1 = %i, N2 = %i, stepsize = %i \n", N1, N2, stepsize);
	// Rprintf("min(N1, N2) = (%i,%i) \n", std::min(N1,N2), std::max(N1,N2));
	int lbd = std::max(std::min(ksize1, ksize2), std::min(N1 - stepsize, N2 - stepsize));
	int ubd = std::max(N1, N2) + stepsize;
	// Rprintf("we consider values of N in [%i,%i] \n", lbd, ubd);
	NumericMatrix w(3, ubd - lbd + 1);
	for (int i = 0; i < w.ncol() ; i++){
	  int m = lbd + i ;
	  w(0,i) = m;
	  if (m >= ksize1){
	    double lcombo = 0;
	    for (int xx = m ; xx > m - ksize1; xx--){
	      lcombo += log(xx);
	    }
	    w(1,i) = lcombo - (n + g) * log(m);
	  }else{
	    w(1,i) = R_NegInf;
	  }
	  if (m >= ksize2){
	    double lcombo = 0;
	    for (int xx = m ; xx > m - ksize2; xx--){
	      lcombo += log(xx);
	    }
	    w(2,i) = lcombo - (n + g) * log(m);
	  }else{
	    w(2,i) = R_NegInf;
	  }
	  // Rprintf("population size N = %i, logweights1  = %f , logweights2 = %f, \n", m, w(1,i), w(2,i));
	}
	// normalize the weights
	double maxlogweights = max(w(1,_));
	double sumweights = 0.0;
	for (int i = 0; i < ubd - lbd + 1; i++){
	  w(1,i) -= maxlogweights;
	  w(1,i) = exp(w(1,i));
	  sumweights += w(1,i);
	}
	for (int i = 0; i < ubd - lbd + 1; i++){
	  w(1,i) = w(1,i) / sumweights;
	}
	maxlogweights = max(w(2,_));
	sumweights = 0.0;
	for (int i = 0; i < ubd - lbd + 1; i++){
	  w(2,i) -= maxlogweights;
	  w(2,i) = exp(w(2,i));
	  sumweights += w(2,i);
	}
	for (int i = 0; i < ubd - lbd + 1; i++){
	  w(2,i) = w(2,i) / sumweights;
	}
	return w;
}