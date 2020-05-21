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
NumericVector compute_psampq_cpp(const int & j,
	const List & clustering,
	const IntegerVector & eta, 
	const NumericMatrix & clusterloglikelihoods,
	const List & theta,
	const IntegerMatrix & V,
	const NumericMatrix & alpha,
	const int & N) {
	int ksize = clustering["ksize"];
	IntegerVector clsize = clustering["clsize"];
	IntegerMatrix clmembers = clustering["clmembers"];
	// the number of clusters excluding record j 
	if (clsize[eta[j]] == 1){
		ksize = ksize - 1 ;
	}
	int n = V.nrow();
	int p = V.ncol();
	NumericMatrix logprod (n,p);
	std::fill(logprod.begin(), logprod.end(), 0.0);
	NumericMatrix lprob_withj (n,p);
	std::fill(lprob_withj.begin(), lprob_withj.end(), NA_REAL);
	NumericMatrix lprob_withoutj (n, p);
	std::fill(lprob_withoutj.begin(), lprob_withoutj.end(), NA_REAL);
	NumericVector lpsampq (n, 0.0);
	NumericVector psampq (n, 1.0);
	for (int l = 0; l < p; l++){
		NumericVector theta_field = theta[l];
		for (int q = 0; q < n; q++){
			// if cluster_q was empty or cluster_q was a singleton cluster with j only
			// just need to compute the products
			if (clsize[q] == 0 || (clsize[q] == 1 && eta[j] == q) ){
				lpsampq[q] += log(theta_field[V(j,l)]);
			}
			// if existing cluster and j was in this cluster and after j is excluded it is still a cluster
			// revert the recursion to get without j probability
			if (clsize[q] > 1 && eta[j] == q){
				lprob_withj(q,l) =  clusterloglikelihoods(q,l);
				lprob_withoutj(q,l) = 0.0;
				logprod(q,l) = 0.0;
				// for all the other members
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember); // index of the member
					if (i != j){// if i is another member in the same cluster with j 
						logprod(q,l) += log( (1 - alpha(q,l)) * (V(i,l) == V(j,l)) + alpha(q,l) * theta_field[V(i,l)]);
						// Rprintf("cluster %i, field %i, imember %i, i %i \n", q, l, imember, i);
					}
				}
				// Rprintf("cluster %i, field %i, logprod %.2f, log with j %.2f \n", q, l, logprod(q,l), lprob_withj(q,l));
				// the unkown is without 
				// pwith = a*p*pwithout + (1-a)*p*prod
				// pwithout = (pwith - (1-a) * p * prod ) / (a*p)
				// psamp = pwith / pwithout 
				// psamp = (a*p*pwithout + (1-a)*p*prod) / pwithout  
				// psamp = pwith * (a * p ) /  (pwith - (1-a) * p * prod )
				logprod(q,l) += log(theta_field[V(j,l)]);
				logprod(q,l) += log(1 - alpha(q,l));
				// we must have pwith - logprod * (1 - a) * p > 0
				// otherwise either pwith is wrong 
				// or logprod is wrong
				if(logprod(q,l) < lprob_withj(q,l)){
					double max_logs = lprob_withj(q,l);
					lprob_withoutj(q,l) =  + max_logs + log(- exp(logprod(q,l) - max_logs) + exp(lprob_withj(q,l) - max_logs)) - log(alpha(q,l)) - log(theta_field[V(j,l)]); 
				}else{
					Rprintf("will come back to this error \n");
					lprob_withoutj(q,l) = lprob_withj(q,l);
				}
				lpsampq[q] += lprob_withj(q,l) - lprob_withoutj(q,l);
				// Rprintf("cluster %i, field %i, logprod %.2f, log without %.2f \n", q, l, logprod(q,l), lprob_withoutj(q,l));
			}
			// if existing cluster and did not include j before 
			// one more step of the recursion
			if (clsize[q] > 0 && eta[j] != q){
				lprob_withoutj(q,l) = clusterloglikelihoods(q,l);
				lprob_withj(q,l) = 0.0;
				logprod(q,l) = 0.0;
				for (int imember = 0; imember < clsize[q]; imember++){
					int i = clmembers(q, imember);
					// must have i != j because they are not in the same cluster
					logprod(q,l) += log( (1 - alpha(q,l)) * (V(i,l) ==  V(j,l)) + alpha(q,l) * theta_field[V(i,l)]);
				}
				// Rprintf("cluster %i, field %i, logprod %.2f, log without j %.2f \n", q, l, logprod(q,l), lprob_withoutj(q,l));
				// the unkown is pwith 
				// pwith = a*p*pwithout + (1-a)*p*prod
				// pwithout = (pwith - (1-a) * p * prod ) / (a*p)
				// psamp = pwith / pwithout 
				// psamp = (a*p*pwithout + (1-a)*p*prod) / pwithout  
				// psamp = pwith * (a * p ) /  (pwith - (1-a) * p * prod )
				logprod(q,l) += log(theta_field[V(j,l)]);
				logprod(q,l) += log(1 - alpha(q,l));
				double max_logs = std::max(logprod(q,l), lprob_withoutj(q,l) + log(alpha(q,l))+ log(theta_field[V(j,l)]));
				lprob_withj(q,l) = max_logs + log(exp(logprod(q,l) - max_logs) + exp(lprob_withoutj(q,l) + log(alpha(q,l)) + log(theta_field[V(j,l)]) - max_logs));
				lpsampq[q] += lprob_withj(q,l) - lprob_withoutj(q,l);
			}
		}
		// Rprintf("cluster %i has ipsamq = %.5f \n", q, lpsampq[q]);
	}
	for (int q = 0; q < n ; q++){
		if (clsize[q] == 0 || (clsize[q] == 1 && eta[j] == q)){
			lpsampq[q] =  lpsampq[q] +  log(double( N - ksize));
		}
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