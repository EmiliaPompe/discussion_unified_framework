// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// compute_loglikelihood_one_cluster_one_field_cpp
double compute_loglikelihood_one_cluster_one_field_cpp(int l, int icluster, const List& clustering, const NumericVector& theta_l, const NumericVector& logtheta_l, const IntegerMatrix& V, const double a);
RcppExport SEXP _couplingdeduplication_compute_loglikelihood_one_cluster_one_field_cpp(SEXP lSEXP, SEXP iclusterSEXP, SEXP clusteringSEXP, SEXP theta_lSEXP, SEXP logtheta_lSEXP, SEXP VSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type icluster(iclusterSEXP);
    Rcpp::traits::input_parameter< const List& >::type clustering(clusteringSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta_l(theta_lSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logtheta_l(logtheta_lSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglikelihood_one_cluster_one_field_cpp(l, icluster, clustering, theta_l, logtheta_l, V, a));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglikelihood_all_clusters_all_fields_cpp
NumericMatrix compute_loglikelihood_all_clusters_all_fields_cpp(const List& clustering, const List& theta, const List& logtheta, const IntegerMatrix& V, const NumericMatrix& alpha);
RcppExport SEXP _couplingdeduplication_compute_loglikelihood_all_clusters_all_fields_cpp(SEXP clusteringSEXP, SEXP thetaSEXP, SEXP logthetaSEXP, SEXP VSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type clustering(clusteringSEXP);
    Rcpp::traits::input_parameter< const List& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type logtheta(logthetaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglikelihood_all_clusters_all_fields_cpp(clustering, theta, logtheta, V, alpha));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglikelihood_all_clusters_one_field_cpp
NumericVector compute_loglikelihood_all_clusters_one_field_cpp(int l, const List& clustering, const NumericVector& theta_l, const NumericVector& logtheta_l, const IntegerMatrix& V, const NumericMatrix& alpha);
RcppExport SEXP _couplingdeduplication_compute_loglikelihood_all_clusters_one_field_cpp(SEXP lSEXP, SEXP clusteringSEXP, SEXP theta_lSEXP, SEXP logtheta_lSEXP, SEXP VSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< const List& >::type clustering(clusteringSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta_l(theta_lSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logtheta_l(logtheta_lSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglikelihood_all_clusters_one_field_cpp(l, clustering, theta_l, logtheta_l, V, alpha));
    return rcpp_result_gen;
END_RCPP
}
// coupled_update_eta
List coupled_update_eta(IntegerVector& eta1, IntegerVector& eta2, List& clustering1, List& clustering2, const NumericMatrix& previous_clusterloglikelihoods1, const NumericMatrix& previous_clusterloglikelihoods2, const List& theta1, const List& theta2, const List& logtheta1, const List& logtheta2, const IntegerMatrix& V, const NumericMatrix& alpha1, const NumericMatrix& alpha2, int N1, int N2, double updateprobability);
RcppExport SEXP _couplingdeduplication_coupled_update_eta(SEXP eta1SEXP, SEXP eta2SEXP, SEXP clustering1SEXP, SEXP clustering2SEXP, SEXP previous_clusterloglikelihoods1SEXP, SEXP previous_clusterloglikelihoods2SEXP, SEXP theta1SEXP, SEXP theta2SEXP, SEXP logtheta1SEXP, SEXP logtheta2SEXP, SEXP VSEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP N1SEXP, SEXP N2SEXP, SEXP updateprobabilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type eta1(eta1SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type eta2(eta2SEXP);
    Rcpp::traits::input_parameter< List& >::type clustering1(clustering1SEXP);
    Rcpp::traits::input_parameter< List& >::type clustering2(clustering2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type previous_clusterloglikelihoods1(previous_clusterloglikelihoods1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type previous_clusterloglikelihoods2(previous_clusterloglikelihoods2SEXP);
    Rcpp::traits::input_parameter< const List& >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< const List& >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const List& >::type logtheta1(logtheta1SEXP);
    Rcpp::traits::input_parameter< const List& >::type logtheta2(logtheta2SEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< int >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< int >::type N2(N2SEXP);
    Rcpp::traits::input_parameter< double >::type updateprobability(updateprobabilitySEXP);
    rcpp_result_gen = Rcpp::wrap(coupled_update_eta(eta1, eta2, clustering1, clustering2, previous_clusterloglikelihoods1, previous_clusterloglikelihoods2, theta1, theta2, logtheta1, logtheta2, V, alpha1, alpha2, N1, N2, updateprobability));
    return rcpp_result_gen;
END_RCPP
}
// init_clustering_cpp
List init_clustering_cpp(const IntegerVector& eta);
RcppExport SEXP _couplingdeduplication_init_clustering_cpp(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(init_clustering_cpp(eta));
    return rcpp_result_gen;
END_RCPP
}
// multinomial_
int multinomial_(const NumericVector& logw, double uniform);
RcppExport SEXP _couplingdeduplication_multinomial_(SEXP logwSEXP, SEXP uniformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type logw(logwSEXP);
    Rcpp::traits::input_parameter< double >::type uniform(uniformSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomial_(logw, uniform));
    return rcpp_result_gen;
END_RCPP
}
// coupled_multinomial_
IntegerVector coupled_multinomial_(const NumericVector& logw1, const NumericVector& logw2, double uniform1, double uniform2);
RcppExport SEXP _couplingdeduplication_coupled_multinomial_(SEXP logw1SEXP, SEXP logw2SEXP, SEXP uniform1SEXP, SEXP uniform2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type logw1(logw1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type logw2(logw2SEXP);
    Rcpp::traits::input_parameter< double >::type uniform1(uniform1SEXP);
    Rcpp::traits::input_parameter< double >::type uniform2(uniform2SEXP);
    rcpp_result_gen = Rcpp::wrap(coupled_multinomial_(logw1, logw2, uniform1, uniform2));
    return rcpp_result_gen;
END_RCPP
}
// update_eta
List update_eta(IntegerVector& eta, List& clustering, const NumericMatrix previous_clusterloglikelihoods, const List& theta, const List& logtheta, const IntegerMatrix& V, const NumericMatrix& alpha, int N, double updateprobability);
RcppExport SEXP _couplingdeduplication_update_eta(SEXP etaSEXP, SEXP clusteringSEXP, SEXP previous_clusterloglikelihoodsSEXP, SEXP thetaSEXP, SEXP logthetaSEXP, SEXP VSEXP, SEXP alphaSEXP, SEXP NSEXP, SEXP updateprobabilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< List& >::type clustering(clusteringSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type previous_clusterloglikelihoods(previous_clusterloglikelihoodsSEXP);
    Rcpp::traits::input_parameter< const List& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type logtheta(logthetaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type updateprobability(updateprobabilitySEXP);
    rcpp_result_gen = Rcpp::wrap(update_eta(eta, clustering, previous_clusterloglikelihoods, theta, logtheta, V, alpha, N, updateprobability));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_couplingdeduplication_compute_loglikelihood_one_cluster_one_field_cpp", (DL_FUNC) &_couplingdeduplication_compute_loglikelihood_one_cluster_one_field_cpp, 7},
    {"_couplingdeduplication_compute_loglikelihood_all_clusters_all_fields_cpp", (DL_FUNC) &_couplingdeduplication_compute_loglikelihood_all_clusters_all_fields_cpp, 5},
    {"_couplingdeduplication_compute_loglikelihood_all_clusters_one_field_cpp", (DL_FUNC) &_couplingdeduplication_compute_loglikelihood_all_clusters_one_field_cpp, 6},
    {"_couplingdeduplication_coupled_update_eta", (DL_FUNC) &_couplingdeduplication_coupled_update_eta, 16},
    {"_couplingdeduplication_init_clustering_cpp", (DL_FUNC) &_couplingdeduplication_init_clustering_cpp, 1},
    {"_couplingdeduplication_multinomial_", (DL_FUNC) &_couplingdeduplication_multinomial_, 2},
    {"_couplingdeduplication_coupled_multinomial_", (DL_FUNC) &_couplingdeduplication_coupled_multinomial_, 4},
    {"_couplingdeduplication_update_eta", (DL_FUNC) &_couplingdeduplication_update_eta, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_couplingdeduplication(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
