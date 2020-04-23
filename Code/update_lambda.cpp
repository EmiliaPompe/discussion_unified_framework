#include <Rcpp.h>
using namespace Rcpp;

// function to update the lambda vector
// where lambda[ilambda] is a label referring to one of N clusters
// takes as arguments 
// * lambda_previous, the vector of current lambdas
// * clustering_previous, the associated clustering
// * ll_previous, the log-likelihoods associated with each cluster
// * p
// * V
// * dimV
// * a
// * N
// [[Rcpp::export]]
List update_lambda(const IntegerVector lambda_previous,
                   const List  clustering_previous,
                   const NumericMatrix ll_previous,
                   const NumericVector & p,
                   const IntegerMatrix & V,
                   const IntegerVector & dimV,
                   const NumericMatrix & a,
                   int N) {
  IntegerVector clsize = clone(wrap(clustering_previous["clsize"]));
  IntegerMatrix clmembers = clone(wrap(clustering_previous["clmembers"]));
  IntegerVector ksize_vec = clone(wrap(clustering_previous["ksize"]));
  int ksize = ksize_vec(0);
  //
  int n = V.nrow();
  int H = V.ncol();
  NumericMatrix clusterloglikelihoods = clone(wrap(ll_previous));
  // 
  int cumdime[H]; 
  cumdime[0] = 0;
  for (int l=1; l<H; l++){
    cumdime[l] = cumdime[l-1] + dimV[l-1];
  }
  //
  IntegerVector lambda = clone(wrap(lambda_previous));
  // loop over components of lambda
  for (int ilambda = 0; ilambda < n; ilambda++){
    int label = lambda[ilambda];
    // first, remove lambda[ilambda] from current partition
    // i.e. translate -1's towards the left
    bool shift = false;
    for (int icol = 0; icol < clsize[label]; icol++){
      if (clmembers(label,icol) == ilambda){
        shift = true;
      }
      if (shift){
        // move entry of next column to current column
        if (icol + 1 < n){
          clmembers(label,icol) = clmembers(label,icol+1);
        } else {
          clmembers(label,icol) = -1;
        }
      } 
    }
    clsize[label] --;
    // if cluster is now empty
    if (clsize[label] == 0){
      // decrement cluster counter
      ksize --;
      // remove log-likelihoods
      std::fill(clusterloglikelihoods.row(label).begin(), clusterloglikelihoods.row(label).end(), 0.);
    } else {
      // need to update log-likelihood associated with what's left in the non-empty cluster
      // for this we use equation (6), which, rearranged, give the likelihood
      // of a cluster with one less member, based on the likelihood of full cluster (reading the recursion backward)
      // loop over fields
      for (int l = 0; l < H; l++){
        // compute second line of equation (6)
        double logprod = 0.;
        for (int othermember = 0; othermember < clsize[label]; othermember ++){
          int iother = clmembers(label,othermember);
          logprod += log((1 - a(label,l)) * (V(ilambda,l) == V(iother,l)) + a(label,l) * p[cumdime[l] + V(iother,l)]);
        }
        logprod += log((1 - a(label,l)) * p[cumdime[l] + V(ilambda,l)]);
        // next we need to define exp(clusterloglikelihoods(l)) - exp(logprod)
        // double max_logs = std::max(logprod, clusterloglikelihoods(label,l));
        double max_logs = 0.;
        clusterloglikelihoods(label,l) = max_logs + log(exp(clusterloglikelihoods(label,l) - max_logs) - exp(logprod - max_logs));
        clusterloglikelihoods(label,l) -= log(a(label,l) * p[cumdime[l] + V(ilambda,l)]);
      }
    }
    // now we have a clustering and associated log likelihoods as if we did not have lambda[ilambda]
    // next we can compute the likelihood associated with a new block with that label in it
    NumericVector newblock_loglikelihood(H, 0.);
    for (int l = 0; l < H; l++){
      newblock_loglikelihood(l) = log(p[cumdime[l] + V(ilambda,l)]);
    }
    NumericMatrix uponlambdajoining_loglikelihood(n,H);
    // vector to aggregate likelihood and prior probabilities for new label
    NumericVector logproba_lambda(n);
    std::fill(logproba_lambda.begin(), logproba_lambda.end(), 0.);
    // next we loop over clusters
    for (int icluster = 0; icluster < n; icluster ++){
      // first compute likelihoods of joining
      if (clsize[icluster]==0){
        // this would be a new cluster
        uponlambdajoining_loglikelihood.row(icluster) = newblock_loglikelihood;
      } else {
        // this would be an existing cluster
        for (int l = 0; l < H; l++){
          // first part of recursion
          uponlambdajoining_loglikelihood(icluster, l) = log(a(icluster,l) * p[cumdime[l] + V(ilambda,l)]) + 
            clusterloglikelihoods(icluster,l);
          double logprod = 0.;
          // next, implement second part of recursion
          for (int clustermember = 0; clustermember < clsize[icluster]; clustermember ++){
            int qprime = clmembers(icluster, clustermember);
            logprod += log((1 - a(icluster,l)) * (V(ilambda,l) == V(qprime,l)) + a(icluster,l) * p[cumdime[l] + V(qprime,l)]);
          }
          logprod += log((1 - a(icluster,l)) * p[cumdime[l] + V(ilambda,l)]);
          double max_logs = std::max(logprod, uponlambdajoining_loglikelihood(icluster, l));
          uponlambdajoining_loglikelihood(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponlambdajoining_loglikelihood(icluster, l) - max_logs));
        }
      }
      // then aggregate priors and conditional likelihood
      if (clsize[icluster] == 0){
        // i.e. P(lambda[ilambda] = q) where q is the label of a new block
        logproba_lambda[icluster] = sum(uponlambdajoining_loglikelihood.row(icluster));
        logproba_lambda[icluster] += (log(N-ksize) - log(n-ksize)); // (this is the prior)
        // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster 
      } else {
        // i.e. P(lambda[ilambda] = q) where q is the label of an existing cluster
        logproba_lambda[icluster] = sum(uponlambdajoining_loglikelihood.row(icluster)) - sum(clusterloglikelihoods.row(icluster));
        logproba_lambda[icluster] += 0.;
      }
    }
    // sample from categorical/multinomial given logproba_lambda
    double maxlogproba_lambda = Rcpp::max(logproba_lambda);
    NumericVector weights = exp(logproba_lambda - maxlogproba_lambda);
    int draw;
    NumericVector cumsumw = cumsum(weights);
    GetRNGstate();
    NumericVector uniform = runif(1);
    PutRNGstate();
    double u = uniform(0) * cumsumw(n-1);
    int running_index = 0;
    double sumw = cumsumw(running_index);
    if (u < sumw){
      draw = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw(running_index);
        running_index ++;
      }
      draw = running_index-1;
    }
    // ... phew that was a bit painful! 
    // thus we update lambda, and the clustering, and associated loglikelihoods appropriately
    lambda[ilambda] = draw;
    // update cluster log likelihoods
    clusterloglikelihoods.row(draw) = uponlambdajoining_loglikelihood.row(draw);
    // update clustering...
    // change number of cluster if need be
    if (clsize[draw] == 0){ 
      ksize += 1;
    }
    // adds index in matrix storing members of each cluster
    clmembers(draw, clsize[draw]) = ilambda;
    // increment cluster size
    clsize[draw] += 1;
  }
  List clustering = List::create(Named("clsize") = clsize, Named("ksize") = ksize, Named("clmembers") = clmembers);  
  return List::create(Named("lambda") = lambda, 
                      Named("clusterloglikelihoods") = clusterloglikelihoods, 
                      Named("clustering") = clustering);
}

