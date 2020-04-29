#include <Rcpp.h>
using namespace Rcpp;

// function to compute the probabitliy of including j and excluding j in each cluster
// the return will serve as an input for couplings
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p which is "theta" in the paper
// V which is the observations
// dimV which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
List compute_include_exclude(const int  j,
	const List  clustering,
	const IntegerVector  eta, 
	const NumericMatrix  clusterloglikelihoods,
	const NumericVector & theta,
	const IntegerMatrix & V,
	const IntegerVector & dimV,
	const NumericMatrix & alpha,
	const int  N) {

	IntegerVector clsize = clone(wrap(clustering["clsize"]));
  	IntegerMatrix clmembers = clone(wrap(clustering["clmembers"]));
  	IntegerVector ksize_vec = clone(wrap(clustering["ksize"]));
  	int ksize = ksize_vec(0);
	// the number of clusters excluding record j 
	if (clsize[eta[j]] == 1){
		ksize = ksize - 1 ;
	}

	int n = V.nrow();
	int L = V.ncol();
	int cumdime[L]; 
	cumdime[0] = 0;
	for (int l=1; l<L; l++){
		cumdime[l] = cumdime[l-1] + dimV[l-1];
	}  
	
	NumericMatrix logprob_include(n,L);
	NumericMatrix logprob_exclude(n,L);		
	NumericVector logprob_newcluster(n);	
	NumericVector lpsampq (n, 0.0);
	NumericVector psampq (n, 1.0);
	std::fill(logprob_newcluster.begin(), logprob_newcluster.end(), 0.0);

	for (int q = 0; q < n; q++){
		lpsampq[q] = 0.0;
		// if cluster_q was empty or cluster_q was a singleton cluster with j only
		if (clsize[q] == 0 || (clsize[q] == 1 && eta[j] == q) ){
			for (int l = 0; l < L ; l++ ){
				logprob_include(q,l) = log(theta[cumdime[l] + V(j,l)]);
				logprob_exclude(q,l) = 0.0;
				lpsampq[q] += log(theta[cumdime[l] + V(j,l)]);
			}
			logprob_newcluster[q] = log( double( N - ksize))- log( double (n - ksize));
			lpsampq[q] =  lpsampq[q] +  log( double( N - ksize))- log( double (n - ksize));
		}
		// if existing cluster and j was in this cluster and after j is excluded it is still a cluster
		// revert the recursion to get without j probability
		if (clsize[q] > 1 && eta[j] == q){
			for (int l = 0; l < L; l++){
				double lprob_withj =  clusterloglikelihoods(q,l);
				double lprob_withoutj = 0.0;
				double logprod = 0.0; // this is the second part of recursion
				// for all the other members
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember); // index of the member
					if (i != j){// if i is another member in the same cluster with j 
						logprod += log( (1 - alpha(q,l)) * (V(i,l) == V(j,l)) + alpha(q,l) * theta[cumdime[l] + V(i,l)]);
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
				logprod += log(theta[cumdime[l] + V(j,l)]);
				logprod += log( 1 - alpha(q,l));
				// we must have pwith - logprod * (1 - a) * p > 0
				// otherwise either pwith is wrong 
				// or logprod is wrong
				if( logprod <= lprob_withj){
					lprob_withoutj =  log(exp(lprob_withj) - exp(logprod)) - log(alpha(q,l)) - log(theta[cumdime[l] + V(j,l)]); 
				}else{
					Rprintf("will come back to this error \n");
					lprob_withoutj = lprob_withj;
				}
				logprob_include(q,l) = lprob_withj;
				logprob_exclude(q,l) = lprob_withoutj;
				lpsampq[q] += lprob_withj - lprob_withoutj;
			}
		}
		// if existing cluster and did not include j before 
		// one more step of the recursion
		if (clsize[q] > 0 && eta[j] != q){
			for (int l = 0; l < L ; l++){
				double lprob_withoutj = clusterloglikelihoods(q,l);
				double lprob_withj = 0.0;
				double logprod = 0.0;
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember);
					// must have i != j because they are not in the same cluster
					logprod += log( (1 - alpha(q,l)) * (V(i,l) ==  V(j,l)) + alpha(q,l) * theta[cumdime[l] + V(i,l)]);
				}
				// the unkown is pwith 
				// pwith = a*p*pwithout + (1-a)*p*prod
				// pwithout = (pwith - (1-a) * p * prod ) / (a*p)
				// psamp = pwith / pwithout 
				// psamp = (a*p*pwithout + (1-a)*p*prod) / pwithout  
				// psamp = pwith * (a * p ) /  (pwith - (1-a) * p * prod )
				logprod += log(theta[cumdime[l] + V(j,l)]);
				logprod += log( 1 - alpha(q,l));
				double max_logs = std::max(logprod, lprob_withoutj + log(alpha(q,l)) + log(theta[cumdime[l] + V(j,l)]));
				lprob_withj = max_logs + log(exp(logprod - max_logs) + exp(lprob_withoutj + log(alpha(q,l)) + log(theta[cumdime[l] + V(j,l)]) - max_logs));
				logprob_include(q,l) = lprob_withj;
				logprob_exclude(q,l) = lprob_withoutj;
				lpsampq[q] += lprob_withj - lprob_withoutj;
			}
		}
	}
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
	return List::create( Named("logprob_include") = logprob_include, Named("logprob_exclude") = logprob_exclude, Named("logprob_newcluster") = logprob_newcluster, Named("psampq") = psampq);
}