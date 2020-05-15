// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // // function to update the eta vector
// // // where eta[ieta] is a label referring to one of N clusters
// // // takes as arguments 
// // // * eta_previous, the vector of current etas
// // // * clustering_previous, the associated clustering
// // // * ll_previous, the log-likelihoods associated with each cluster
// // // * theta, list of frequencies for each field
// // // * V, observations
// // // * dimV, number of possibilities for each field
// // // * alpha, corruption probabilities for each field and cluster
// // // * N
// // // [[Rcpp::export]]
// // List update_eta_cpp(const IntegerVector eta_previous,
// //                 const List  clustering_previous,
// //                 const NumericMatrix ll_previous,
// //                 const List & theta,
// //                 const IntegerMatrix & V,
// //                 const IntegerVector & dimV,
// //                 const NumericMatrix & alpha,
// //                 int N,
// //                 const NumericVector & beta0, double s) {
// //   IntegerVector clsize = clone(wrap(clustering_previous["clsize"]));
// //   IntegerMatrix clmembers = clone(wrap(clustering_previous["clmembers"]));
// //   IntegerVector ksize_vec = clone(wrap(clustering_previous["ksize"]));
// //   NumericMatrix new_alpha = clone(wrap(alpha));
// //   int ksize = ksize_vec(0);
// //   //
// //   int n = V.nrow();
// //   int p = V.ncol();
// //   NumericMatrix clusterloglikelihoods = clone(wrap(ll_previous));
// //   // 
// //   IntegerVector eta = clone(wrap(eta_previous));
// //   // loop over components of eta
// //   for (int ieta = 0; ieta < n; ieta++){
// //     int label = eta[ieta];
// //     // first, remove eta[ieta] from current partition
// //     // i.e. translate -1's towards the left
// //     bool shift = false;
// //     for (int icol = 0; icol < clsize[label]; icol++){
// //       if (clmembers(label,icol) == ieta){
// //         shift = true;
// //       }
// //       if (shift){
// //         // move entry of next column to current column
// //         if (icol + 1 < n){
// //           clmembers(label,icol) = clmembers(label,icol+1);
// //         } else {
// //           clmembers(label,icol) = -1;
// //         }
// //       } 
// //     }
// //     clsize[label] --;
// //     // if cluster is now empty
// //     if (clsize[label] == 0){
// //       // decrement cluster counter
// //       ksize --;
// //       // remove log-likelihoods associated with that cluster
// //       std::fill(clusterloglikelihoods.row(label).begin(), clusterloglikelihoods.row(label).end(), NA_REAL);
// //       // remove corresponding alpha
// //       // std::fill(new_alpha.row(label).begin(), new_alpha.row(label).end(), R_NegInf);
// //     } else {
// //       // otherwise need to update log-likelihood associated with what's left in the non-empty cluster
// //       // for this we use equation (6), which, rearranged, gives the likelihood
// //       // of a cluster with one less member, based on the likelihood of full cluster (reading the recursion backward)
// //       // loop over fields
// //       // for (int l = 0; l < p; l++){
// //       //   double eval_ = 0.;
// //       //   eval_ = clusterloglikelihoods(label,l);
// //       //   NumericVector theta_field = theta[l];
// //       //   // compute second line of equation (6)
// //       //   double logprod = 0.;
// //       //   for (int othermember = 0; othermember < clsize[label]; othermember ++){
// //       //     int iother = clmembers(label,othermember);
// //       //     logprod += log((1 - alpha(label,l)) * (V(ieta,l) == V(iother,l)) + alpha(label,l) * theta_field(V(iother,l)));
// //       //   }
// //       //   logprod += log((1 - alpha(label,l)) * theta_field(V(ieta,l)));
// //       //   // next we need to define exp(clusterloglikelihoods(l)) - exp(logprod)
// //       //   double max_logs = std::max(logprod, eval_);
// //       //   double eval_2 = max_logs + log(exp(eval_ - max_logs) - exp(logprod - max_logs));
// //       //   eval_2 -= log(alpha(label,l) * theta_field(V(ieta,l)));
// //       //   if (traits::is_nan<REALSXP>(eval_2)){
// //       //     // clusterloglikelihoods(label,l) = R_NegInf;
// //       //     // in that case we recompute 
// //       //   } else {
// //       //     clusterloglikelihoods(label,l) = eval_2;
// //       //   }
// //       // }
// //       // alternatively we can recompute likelihood of cluster from scratch,
// //       // which seems more numerically stable
// //       NumericVector cl_likelihood_field(p);
// //       std::fill(cl_likelihood_field.begin(), cl_likelihood_field.end(), 0.0);
// //       // compute likelihood recursively over members of cluster
// //       // first member (associated row in V)
// //       int j = clmembers(label,0);
// //       for (int l = 0; l < p; l++){
// //         NumericVector theta_field = theta[l];
// //         // recall p[cumdime[l] + V(j,l)] is theta_{l, v(j,l)} in the paper
// //         cl_likelihood_field(l) += log(theta_field(V(j,l)));
// //       }
// //       if (clsize[label] > 1){
// //         // if more members, implement recursion
// //         // loop over remaining members
// //         for (int imember = 1; imember < clsize[label]; imember++){
// //           // original index of that member (i.e. corresponding row in V)
// //           int q = clmembers(label,imember);
// //           // loop over fields 
// //           for (int l = 0; l < p; l++){
// //             NumericVector theta_field = theta[l];
// //             double logprod = 0.;
// //             // to implement first part of recursion
// //             // multiply by associated alpha' * V(q,l)
// //             cl_likelihood_field(l) += log(new_alpha(label,l) * theta_field(V(q,l)));
// //             // next, implement second part of recursion
// //             for (int othermember = 0; othermember < imember; othermember ++){
// //               int qprime = clmembers(label,othermember);
// //               int sameV = 0;
// //               if (V(q,l) == V(qprime,l)){
// //                 sameV = 1;
// //               }
// //               logprod += log((1 - new_alpha(label,l)) * sameV + new_alpha(label,l) * theta_field(V(qprime,l)));
// //             }
// //             logprod += log((1 - new_alpha(label,l)) * theta_field(V(q,l)));
// //             // next we need to define exp(logprod) + exp(cl_likelihood_field(l))
// //             double max_logs = std::max(logprod, cl_likelihood_field(l));
// //             cl_likelihood_field(l) = max_logs + log(exp(logprod - max_logs) + exp(cl_likelihood_field(l) - max_logs));
// //           }
// //         }
// //       }
// //       for (int l = 0; l < p; l++){
// //         clusterloglikelihoods(label,l) = cl_likelihood_field(l);
// //       }
// //     }
// //     // now we have a clustering and associated log likelihoods as if we did not have eta[ieta] in the partition
// //     // next we can compute the likelihood associated with a new block with that label in it
// //     NumericVector newblock_loglikelihood(p, 0.);
// //     for (int l = 0; l < p; l++){
// //       NumericVector theta_field = theta[l];
// //       newblock_loglikelihood(l) = log(theta_field(V(ieta,l)));
// //     }
// //     NumericMatrix uponetajoining_loglikelihood(n,p);
// //     // vector to aggregate likelihood and prior probabilities for new label
// //     NumericVector logproba_eta(n);
// //     std::fill(logproba_eta.begin(), logproba_eta.end(), 0.);
// //     // next we loop over clusters
// //     for (int icluster = 0; icluster < n; icluster ++){
// //       // first compute likelihoods of joining
// //       if (clsize[icluster]==0){
// //         // this would be a new cluster
// //         uponetajoining_loglikelihood.row(icluster) = newblock_loglikelihood;
// //       } else {
// //         // this would be an existing cluster
// //         for (int l = 0; l < p; l++){
// //           NumericVector theta_field = theta[l];
// //           // first part of recursion
// //           uponetajoining_loglikelihood(icluster, l) = log(new_alpha(icluster,l) * theta_field(V(ieta,l))) + 
// //             clusterloglikelihoods(icluster,l);
// //           double logprod = 0.;
// //           // next, implement second part of recursion
// //           for (int clustermember = 0; clustermember < clsize[icluster]; clustermember ++){
// //             int qprime = clmembers(icluster, clustermember);
// //             logprod += log((1 - new_alpha(icluster,l)) * (V(ieta,l) == V(qprime,l)) + new_alpha(icluster,l) * theta_field(V(qprime,l)));
// //           }
// //           logprod += log((1 - new_alpha(icluster,l)) * theta_field(V(ieta,l)));
// //           double max_logs = std::max(logprod, uponetajoining_loglikelihood(icluster, l));
// //           uponetajoining_loglikelihood(icluster, l) = max_logs + log(exp(logprod - max_logs) + exp(uponetajoining_loglikelihood(icluster, l) - max_logs));
// //         }
// //       }
// //       // then aggregate priors and conditional likelihood
// //       if (clsize[icluster] == 0){
// //         // i.e. P(eta[ieta] = q) where q is the label of a new block
// //         logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster));
// //         logproba_eta[icluster] += (log(N-ksize)); // CAREFUL CHECK - log(n-ksize)); // (this is the prior)
// //         // note that when N = ksize this should be -Inf so it becomes impossible to create a new cluster 
// //       } else {
// //         logproba_eta[icluster] = sum(uponetajoining_loglikelihood.row(icluster) - clusterloglikelihoods.row(icluster));
// //         logproba_eta[icluster] += 0.;
// //         // if (traits::is_infinite<REALSXP>(logproba_eta[icluster])){
// //         //   logproba_eta[icluster] = -100000000.0;
// //         // }
// //         // Rcout << "cluster logprob " << logproba_eta[icluster] << "\n";
// //       }
// //     }
// //     
// //     // sample from categorical/multinomial given logproba_eta
// //     double maxlogproba_eta = Rcpp::max(logproba_eta);
// //     // Rcout << maxlogproba_eta << std::endl;
// //     NumericVector weights = exp(logproba_eta - maxlogproba_eta);
// //     int draw = 0;
// //     NumericVector cumsumw = cumsum(weights);
// //     GetRNGstate();
// //     NumericVector uniform = runif(1);
// //     PutRNGstate();
// //     double u = uniform(0) * cumsumw(n-1);
// //     int running_index = 0;
// //     double sumw = cumsumw(running_index);
// //     if (u <= sumw){
// //       draw = running_index;
// //     } else {
// //       while (u > sumw){
// //         sumw = cumsumw(running_index);
// //         running_index ++;
// //       }
// //       draw = running_index-1;
// //     }
// //     // thus we update eta, and the clustering, and associated loglikelihoods appropriately
// //     eta[ieta] = draw;
// //     // update cluster log likelihoods
// //     clusterloglikelihoods.row(draw) = uponetajoining_loglikelihood.row(draw);
// //     // if cluster is new, we need to draw a vector alpha 
// //     // and we need to update the partition
// //     // change number of cluster if need be
// //     if (clsize[draw] == 0){ 
// //       ksize += 1;
// //       // if new cluster, draw new alpha from prior 
// //       // NumericVector zz = rnorm(p);
// //       // zz = beta0 + s * zz;
// //       // zz = exp(zz)/(1+exp(zz));
// //       // new_alpha.row(draw) = zz;
// //     }
// //     // adds index in matrix storing members of each cluster
// //     clmembers(draw, clsize[draw]) = ieta;
// //     // increment cluster size
// //     clsize[draw] += 1;
// //   }
// //   List clustering = List::create(Named("clsize") = clsize, Named("ksize") = ksize, Named("clmembers") = clmembers);  
// //   return List::create(Named("eta") = eta, 
// //                       Named("alpha") = new_alpha, 
// //                       Named("clusterloglikelihoods") = clusterloglikelihoods, 
// //                       Named("clustering") = clustering);
// // }
// 
