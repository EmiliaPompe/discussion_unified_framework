#include <Rcpp.h>
using namespace Rcpp;

// function that actually does the recursion
double ll_cluster_field(int clustersize, 
                        const IntegerVector & clmembers,
                        const NumericVector & theta_l, 
                        const IntegerVector & Vfield, const double a){
  double log_likelihood = NA_REAL ; 
  double logprod,max_logs;
  int q,j,qprime,sameV;
  double loga = log(a);
  double log1minusa = log(1-a);
  if (clustersize == 0){
    // do nothing and return NA
  } else {
    // compute likelihood recursively over members of cluster
    log_likelihood = 0.0;
    // first member (associated row in V)
    j = clmembers(0);
    // recall theta_l[ V(j,field)] is theta_{l, v(j,l)} in the paper
    log_likelihood += log(theta_l[Vfield(j)]);
    if (log_likelihood > 1){
      // if more members, implement recursion
      // loop over remaining members
      for (int imember = 1; imember < clustersize; imember++){
        // original index of that member (i.e. corresponding row in V)
        q = clmembers(imember);
        logprod = 0.;
        // to implement first part of recursion
        // multiply by associated alpha' * theta_l[V(q,l)]
        log_likelihood += loga + log(theta_l[Vfield(q)]);
        // next, implement second part of recursion
        for (int othermember = 0; othermember < imember; othermember ++){
          qprime = clmembers(othermember);
          sameV = 0;
          if (Vfield(q) == Vfield(qprime)){
            sameV = 1;
          }
          logprod += log((1 - a) * sameV + a * theta_l[Vfield(qprime)]);
        }
        logprod += log1minusa + log(theta_l[Vfield(q)]);
        max_logs = std::max(logprod, log_likelihood);
        log_likelihood = max_logs + log(exp(logprod - max_logs) + exp(log_likelihood - max_logs));
      }
    }
  }
  return log_likelihood;  
}

// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times H = ncol(V)
// if a cluster is empty, the corresponding row is filled with zeros
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to 
// values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// a is alpha_ prime[icluster,l] in the paper, and it is alpha[icluster,l] of
// the R code

// [[Rcpp::export]]
double compute_loglikelihood_one_cluster_one_field_cpp(int l ,
                                                       int icluster, 
                                                       const List & clustering,
                                                       const NumericVector & theta_l,
                                                       const IntegerMatrix & V,
                                                       const double a) {
  const IntegerVector & clsize = clustering["clsize"];
  const IntegerMatrix & clmembers = clustering["clmembers"];
  double log_likelihood = NA_REAL ; 
  log_likelihood = ll_cluster_field(clsize[icluster], 
                                    clmembers(icluster,_),
                                    theta_l, 
                                    V(_,l), a);
  return log_likelihood;
}


// function to compute log likelihoods associated with each cluster and each field 
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with NA
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// Mvec which contains the number of possible categories in each column of V
// alpha is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericMatrix compute_loglikelihood_all_clusters_all_fields_cpp(const List & clustering,
                                                                const List & theta,
                                                                const IntegerMatrix & V,
                                                                const NumericMatrix & alpha) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  int p = V.ncol();
  // 
  NumericMatrix clusterloglikelihoods(n, p);
  std::fill(clusterloglikelihoods.begin(), clusterloglikelihoods.end(), NA_REAL);
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      for (int field = 0; field < p; field++){
        clusterloglikelihoods(icluster, field) = ll_cluster_field(clsize[icluster],  clmembers(icluster,_),
                              theta[field], V(_,field), alpha(icluster, field));
      }
    }
  }
  return clusterloglikelihoods;
}

// function to compute log likelihoods associated with each cluster and one particular field
// returns a matrix of size n = nrow(V) times p = ncol(V)
// if a cluster is empty, the corresponding row is filled with NA
// takes as argument 'clustering', as obtained e.g. with "init_clustering"
// p is "theta" in the paper
// V is the observations (careful, values need to start at zero and go to values in Mvec-1)
// mVec which contains the number of possible categories in each column of V
// a is alpha prime in the paper and given in a matrix with n row, one for each cluster, and ncol(V) columns

// [[Rcpp::export]]
NumericVector compute_loglikelihood_all_clusters_one_field_cpp(int l ,
                                                               const List & clustering,
                                                               const NumericVector & theta_l,
                                                               const IntegerMatrix & V,
                                                               const NumericMatrix & alpha) {
  
  // int ksize = clustering["ksize"];
  IntegerVector clsize = clustering["clsize"];
  IntegerMatrix clmembers = clustering["clmembers"];
  int n = V.nrow();
  // 
  NumericVector cl_likelihood_field(n);
  std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), NA_REAL);
  
  // loop over non-empty clusters and compute associated likelihoods
  for (int icluster = 0; icluster < n; icluster ++){
    // if empty cluster, do nothing
    if (clsize[icluster] == 0){
      // do nothing
    } else {
      cl_likelihood_field(icluster) = ll_cluster_field(clsize[icluster],  clmembers(icluster,_),
                          theta_l, V(_,l), alpha(icluster, l));
    }
  }
  return cl_likelihood_field;
}


// function to update the eta vector
// where eta[ieta] is a label referring to one of N clusters
// takes as arguments 
// * eta_previous, the vector of current etas
// * clustering_previous, the associated clustering
// * ll_previous, the log-likelihoods associated with each cluster
// * theta, list of frequencies for each field
// * V, observations
// * dimV, number of possibilities for each field
// * alpha, corruption probabilities for each field and cluster
// * N
// [[Rcpp::export]]
List update_eta_cpp(const IntegerVector eta_previous,
                    const List  clustering_previous,
                    const NumericMatrix ll_previous,
                    const List & theta,
                    const IntegerMatrix & V,
                    const IntegerVector & dimV,
                    const NumericMatrix & alpha,
                    int N,
                    const NumericVector & beta0, double s) {
  IntegerVector clsize = clone(wrap(clustering_previous["clsize"]));
  IntegerMatrix clmembers = clone(wrap(clustering_previous["clmembers"]));
  IntegerVector ksize_vec = clone(wrap(clustering_previous["ksize"]));
  NumericMatrix new_alpha = clone(wrap(alpha));
  int ksize = ksize_vec(0);
  //
  int n = V.nrow();
  int p = V.ncol();
  NumericMatrix clusterloglikelihoods = clone(wrap(ll_previous));
  // 
  IntegerVector eta = clone(wrap(eta_previous));
  // loop over components of eta
  for (int ieta = 0; ieta < n; ieta++){
    int label = eta[ieta];
    // first, remove eta[ieta] from current partition
    // i.e. translate -1's towards the left
    bool shift = false;
    for (int icol = 0; icol < clsize[label]; icol++){
      if (clmembers(label,icol) == ieta){
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
      // remove log-likelihoods associated with that cluster
      std::fill(clusterloglikelihoods.row(label).begin(), clusterloglikelihoods.row(label).end(), NA_REAL);
      // remove corresponding alpha
      // std::fill(new_alpha.row(label).begin(), new_alpha.row(label).end(), R_NegInf);
    } else {
      // recompute likelihood of cluster from scratch,
      // NumericVector cl_likelihood_field(p);
      // std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), 0.0);
      for (int field = 0; field < p; field ++){
        clusterloglikelihoods(label,field) = compute_loglikelihood_one_cluster_one_field_cpp(field,
                            label, clustering_previous, theta[field], V, new_alpha(label,field));
      }        
      
      // // compute likelihood recursively over members of cluster
      // // first member (associated row in V)
      // int j = clmembers(label,0);
      // for (int l = 0; l < p; l++){
      //   NumericVector theta_field = theta[l];
      //   // recall p[cumdime[l] + V(j,l)] is theta_{l, v(j,l)} in the paper
      //   cl_likelihood_field(l) += log(theta_field(V(j,l)));
      // }
      // if (clsize[label] > 1){
      //   // if more members, implement recursion
      //   // loop over remaining members
      //   for (int imember = 1; imember < clsize[label]; imember++){
      //     // original index of that member (i.e. corresponding row in V)
      //     int q = clmembers(label,imember);
      //     // loop over fields 
      //     for (int l = 0; l < p; l++){
      //       NumericVector theta_field = theta[l];
      //       double logprod = 0.;
      //       // to implement first part of recursion
      //       // multiply by associated alpha' * V(q,l)
      //       cl_likelihood_field(l) += log(new_alpha(label,l) * theta_field(V(q,l)));
      //       // next, implement second part of recursion
      //       for (int othermember = 0; othermember < imember; othermember ++){
      //         int qprime = clmembers(label,othermember);
      //         int sameV = 0;
      //         if (V(q,l) == V(qprime,l)){
      //           sameV = 1;
      //         }
      //         logprod += log((1 - new_alpha(label,l)) * sameV + new_alpha(label,l) * theta_field(V(qprime,l)));
      //       }
      //       logprod += log((1 - new_alpha(label,l)) * theta_field(V(q,l)));
      //       // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
      //       double max_logs = std::max(logprod, cl_likelihood_field(l));
      //       cl_likelihood_field(l) = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field(l) - max_logs));
      //     }
      //   }
      // }
      // for (int l = 0; l < p; l++){
      //   clusterloglikelihoods(label,l) = cl_likelihood_field(l);
      // }
    }
    // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
    // next we can compute the likelihood associated with a new block with that label in it
    NumericVector newblock_loglikelihood(p, 0.);
    for (int l = 0; l < p; l++){
      NumericVector theta_field = theta[l];
      newblock_loglikelihood(l) = log(theta_field(V(ieta,l)));
    }
    NumericMatrix uponetajoining_loglikelihood(n,p);
    // vector to aggregate likelihood and prior probabilities for new label
    NumericVector logproba_eta(n);
    std::fill(logproba_eta.begin(), logproba_eta.end(), 0.);
    // next we loop over clusters
    for (int icluster = 0; icluster < n; icluster ++){
      // first compute likelihoods of joining
      if (clsize[icluster]==0){
        // this would be a new cluster
        uponetajoining_loglikelihood.row(icluster) = newblock_loglikelihood;
      } else {
        // this would be an existing cluster
        for (int l = 0; l < p; l++){
          NumericVector theta_field = theta[l];
          // first part of recursion
          uponetajoining_loglikelihood(icluster, l) = log(new_alpha(icluster,l) * theta_field(V(ieta,l))) + 
            clusterloglikelihoods(icluster,l);
          double logprod = 0.;
          // next, implement second part of recursion
          for (int clustermember = 0; clustermember < clsize[icluster]; clustermember ++){
            int qprime = clmembers(icluster, clustermember);
            logprod += log((1 - new_alpha(icluster,l)) * (V(ieta,l) == V(qprime,l)) + new_alpha(icluster,l) * theta_field(V(qprime,l)));
          }
          logprod += log((1 - new_alpha(icluster,l)) * theta_field(V(ieta,l)));
          double max_logs = std::max(logprod, uponetajoining_loglikelihood(icluster, l));
          uponetajoining_loglikelihood(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood(icluster, l) - max_logs));
        }
      }
      // then aggregate priors and conditional likelihood
      if (clsize[icluster] == 0){
        // i.e. P(eta[ieta] = q) where q is the label of a new block
        logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster));
        logproba_eta[icluster] += (log(N-ksize)); // CAREFUL CHECK - log(n-ksize));
        // logproba_eta[icluster] += (log(N-ksize) - log(n-ksize)); // CAREFUL CHECK
        // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster 
      } else {
        logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster) - clusterloglikelihoods.row(icluster));
        logproba_eta[icluster] += 0.;
      }
    }
    
    // sample from categorical/multinomial given logproba_eta
    double maxlogproba_eta = Rcpp::max(logproba_eta);
    // Rcout << maxlogproba_eta << std::endl;
    NumericVector weights = exp(logproba_eta - maxlogproba_eta);
    int draw = 0;
    NumericVector cumsumw = cumsum(weights);
    GetRNGstate();
    NumericVector uniform = runif(1);
    PutRNGstate();
    double u = uniform(0) * cumsumw(n-1);
    int running_index = 0;
    double sumw = cumsumw(running_index);
    if (u <= sumw){
      draw = running_index;
    } else {
      while (u > sumw){
        sumw = cumsumw(running_index);
        running_index ++;
      }
      draw = running_index-1;
    }
    // thus we update eta, and the clustering, and associated loglikelihoods appropriately
    eta[ieta] = draw;
    // update cluster log likelihoods
    clusterloglikelihoods.row(draw) = uponetajoining_loglikelihood.row(draw);
    // if cluster is new, we need to draw a vector alpha 
    // and we need to update the partition
    // change number of cluster if need be
    if (clsize[draw] == 0){ 
      ksize += 1;
      // if new cluster, draw new alpha from prior 
      // NumericVector zz = rnorm(p);
      // zz = beta0 + s * zz;
      // zz = exp(zz)/(1+exp(zz));
      // new_alpha.row(draw) = zz;
    }
    // adds index in matrix storing members of each cluster
    clmembers(draw, clsize[draw]) = ieta;
    // increment cluster size
    clsize[draw] += 1;
  }
  List clustering = List::create(Named("clsize") = clsize, Named("ksize") = ksize, Named("clmembers") = clmembers);  
  return List::create(Named("eta") = eta, 
                      Named("alpha") = new_alpha, 
                      Named("clusterloglikelihoods") = clusterloglikelihoods, 
                      Named("clustering") = clustering);
}

