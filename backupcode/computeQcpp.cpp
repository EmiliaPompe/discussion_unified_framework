#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeQcpp(int j,
	IntegerVector lambda, 
	NumericVector p,
	IntegerMatrix V,
	IntegerVector cumdime,
	NumericMatrix a,
	int n,
	int H,
	int N) {
	// initialize pclusterM, probM and psampq;
	NumericVector psampq (n,1.0);
	IntegerVector clsize (n, 0);
	int ksize = 0;
	NumericVector probM (H, 1.0);
	NumericVector pclusterM (H, 1.0);
	NumericVector cumprod (H, 1.0);

	// compute clsize and ksize 
	for (int i = 0; i < n ; i++){
	  // Rprintf("the %i-th record belongs to cluster %i \n", i, lambda[i]);
		clsize[ lambda[ i ] ] += 1;
	}
	for (int l = 0; l < n;  l++){
		if (clsize[ l ] > 0 ) ksize += 1;
	}
	if (clsize[j] == 1){
		ksize -= 1;
	}
	// for(int i=0; i < clsize.length(); ++i){
	//   Rprintf("the value of v[%i] : %i \n", i, clsize[i]);
	// }
	// Rprintf("there are %i clusters \n", ksize);
	for (int q = 0; q < n; q++){
		// first see if q is a new cluster
		if (clsize[q] == 0){
			for (int l = 0; l < H; l++){
				psampq[q] *= p[cumdime[l]+ V(j,l)]; 
			}
			psampq[q] = psampq[q] * (N - ksize) / (n - ksize);
		}
		else{
			// if q is an observed cluster 
			// psampq is the ratio of include j vs exclude j
			for (int l = 0 ; l < H; l ++){
	    		probM[ l ] = 1.0;
	    		pclusterM[ l ] = 1.0;
	    		cumprod[ l ] = 1.0;
	  		}
	  		// pclusterM is the probability of excluding j in q
	  		// probM is the probability of including j in q
	  		bool firstElement = true;
	  		for (int i = 0; i < n ; i++){
	  			if (lambda[i] == q && j != i){
	  				if (firstElement){
	  					firstElement = false;
	  					for (int l = 0; l < H; l++){
	  						pclusterM[l] = p[cumdime[l] + V(j,l)];
	  					}
	  				}else{
	  					for(int k = 0; k < i; k++){
	  						// for all the inspected record 
	  						if (lambda[k] == q && k != j){
	  							for(int l = 0 ; l < H; l++ ){
	  								cumprod[l] *= (1 - a(q,l)) * (V(k,l)== V(i,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
	  								pclusterM[l] = pclusterM[l] * a(q,l) * p[cumdime[l] + V(i,l)]+ (1 - a(q,l)) * p[cumdime[l] + V(i,l)] * cumprod[l];	
	  							}
	  						}
	  					}
	  				}
	  				for (int l = 0; l < H; l++){
	  					probM[l] *= (1 - a(q,l)) * (V(i,l) == V(j,l)) + a(q,l) * p[cumdime[l] + V(i,l)];
	  				}
	  			}
	  		}
	  		for (int l = 0; l < H; l++){
	  			probM[l] = p[cumdime[l] + V(j,l)] * ( a(q,l) + (1 - a(q,l)) * probM[l] / pclusterM[l]);
	  		}
	  		psampq[q] = 1.0;
	  		for (int l = 0; l < H; l++){
	  			psampq[q] *= probM[l];
	  		}
		}
	}

	// normalize psampq
	double sumpsampq = 0; 
	for (int q = 0; q < n ; q++){
	  sumpsampq += psampq[q];
	}
	for (int q = 0; q < n ; q++){
	  psampq[q] = psampq[q] / sumpsampq;
	}
	return psampq;
}