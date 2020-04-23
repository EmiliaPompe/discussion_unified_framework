#include <Rcpp.h>
using namespace Rcpp;

// function to compute the full conditional distribution of labels of record j
// returns a vector of length n 
// the return will serve as an input for couplings
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p which is "theta" in the paper
// V which is the observations
// dimV which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_psampq(const int & j,
	const List & clustering,
	const IntegerVector & lambda, 
	const NumericMatrix & clusterloglikelihoods,
	const NumericVector & p,
	const IntegerMatrix & V,
	const IntegerVector & dimV,
	const NumericMatrix & a,
	const int & N) {
	int ksize = clustering["ksize"];
	IntegerVector clsize = clustering["clsize"];
	IntegerMatrix clmembers = clustering["clmembers"];
	// the number of clusters excluding record j 
	if (clsize[lambda[j]] == 1){
		ksize = ksize - 1 ;
	}

	int n = V.nrow();
	int H = V.ncol();
	// 
	int cumdime[H]; 
	cumdime[0] = 0;
	for (int l=1; l<H; l++){
		cumdime[l] = cumdime[l-1] + dimV[l-1];
	}  
	NumericVector lpsampq (n, 0.0);
	NumericVector psampq (n, 1.0);
	for (int q = 0; q < n; q++){
		lpsampq[q] = 0.0;
		// if cluster_q was empty or cluster_q was a singleton cluster with j only
		// just need to compute the products
		if (clsize[q] == 0 || (clsize[q] == 1 && lambda[j] == q) ){
			for (int l = 0; l < H ; l++ ){
				lpsampq[q] += log(p[cumdime[l] + V(j,l)]);
			}
			lpsampq[q] =  lpsampq[q] +  log( double( N - ksize))- log( double (n - ksize));
		}
		// if existing cluster and j was in this cluster and after j is excluded it is still a cluster
		// revert the recursion to get without j probability
		if (clsize[q] > 1 && lambda[j] == q){
			for (int l = 0; l < H; l++){
				double lprob_withj =  clusterloglikelihoods(q,l);
				double lprob_withoutj = 0.0;
				double logprod = 0.0; // this is the second part of recursion
				// for all the other members
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember); // index of the member
					if (i != j){// if i is another member in the same cluster with j 
						logprod += log( (1 - a(q,l)) * (V(i,l) == V(j,l)) + a(q,l) * p[cumdime[l] + V(i,l)]);
						// Rprintf("cluster %i, field %i, imember %i, i %i \n", q, l, imember, i);
					}

				}
				// Rprintf("cluster %i, field %i, logprod %.2f, log with j %.2f \n", q, l, logprod, lprob_withj);
				// the unkown is without 
				// pwith = a*p*pwithout + (1-a)*p*prod
				// pwithout = (pwith - (1-a) * p * prod ) / (a*p)
				// psamp = pwith / pwithout 
				// psamp = (a*p*pwithout + (1-a)*p*prod) / pwithout  
				// psamp = pwith * (a * p ) /  (pwith - (1-a) * p * prod )
				logprod += log(p[cumdime[l] + V(j,l)]);
				logprod += log( 1 - a(q,l));
				// we must have pwith - logprod * (1 - a) * p > 0
				// otherwise either pwith is wrong 
				// or logprod is wrong
				if( logprod < lprob_withj){
					double max_logs = lprob_withj;
					lprob_withoutj =  + max_logs + log(- exp(logprod - max_logs) + exp(lprob_withj - max_logs)); 
				}else{
					// Rprintf("will come back to this error \n");
					lprob_withoutj = lprob_withj;
				}
				lpsampq[q] += lprob_withj - lprob_withoutj;
				// Rprintf("cluster %i, field %i, logprod %.2f, log without %.2f \n", q, l, logprod, lprob_withoutj);
			}
		}
		// if existing cluster and did not include j before 
		// one more step of the recursion
		if (clsize[q] > 0 && lambda[j] != q){
			for (int l = 0; l < H ; l++){
				double lprob_withoutj = clusterloglikelihoods(q,l);
				double lprob_withj = 0.0;
				double logprod = 0.0;
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember);
					// must have i != j because they are not in the same cluster
					logprod += log( (1 - a(q,l)) * (V(i,l) ==  V(j,l)) + a(q,l) * p[cumdime[l] + V(i,l)]);
				}
				// Rprintf("cluster %i, field %i, logprod %.2f, log without j %.2f \n", q, l, logprod, lprob_withoutj);
				// the unkown is pwith 
				// pwith = a*p*pwithout + (1-a)*p*prod
				// pwithout = (pwith - (1-a) * p * prod ) / (a*p)
				// psamp = pwith / pwithout 
				// psamp = (a*p*pwithout + (1-a)*p*prod) / pwithout  
				// psamp = pwith * (a * p ) /  (pwith - (1-a) * p * prod )
				logprod += log(p[cumdime[l] + V(j,l)]);
				logprod += log( 1 - a(q,l));
				double max_logs = std::max(logprod, lprob_withoutj + log(a(q,l)) + log(p[cumdime[l] + V(j,l)]));
				lprob_withj = max_logs + log(exp(logprod - max_logs) + exp(lprob_withoutj + log(a(q,l)) + log(p[cumdime[l] + V(j,l)]) - max_logs));
				lpsampq[q] += lprob_withj - lprob_withoutj;
			}
		}
		// Rprintf("cluster %i has ipsamq = %.5f \n", q, lpsampq[q]);
	}
	// normalize the lpsamp to get psampq 
	// double max_logs = max(lpsampq);
	// lpsampq = lpsampq - max_logs;
	// psampq = exp(lpsampq) / sum(exp(lpsampq));
	// for (int q= 0; q < n ; q++ ){
	// 	Rprintf("cluster %i has psamq proportional to %.5f \n", q, exp(lpsampq[q]));
	// }
	double max_logs = max(lpsampq);
	double sumexplogsampq = 0.0;
	for (int q = 0; q < n ; q++){
		lpsampq[q] = lpsampq[q] - max_logs;
		// Rprintf("cluster %i has psamq proportional to %.5f \n", q, exp(lpsampq[q]));
		sumexplogsampq += exp(lpsampq[q]);
	}
	for (int q = 0; q < n ; q++){
		psampq[q] = exp(lpsampq[q]) / sumexplogsampq;
	}
	return psampq;
}