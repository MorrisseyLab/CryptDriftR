// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Biased_Cont_Labelling_Theory_c
Rcpp::List Biased_Cont_Labelling_Theory_c(double kappa, double alpha, double beta, int Ns);
RcppExport SEXP _CryptDriftR_Biased_Cont_Labelling_Theory_c(SEXP kappaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    rcpp_result_gen = Rcpp::wrap(Biased_Cont_Labelling_Theory_c(kappa, alpha, beta, Ns));
    return rcpp_result_gen;
END_RCPP
}
// LinearGillespie
arma::umat LinearGillespie(int numSim, arma::umat nu, arma::mat a_mat, arma::uvec x0, arma::vec time_vec);
RcppExport SEXP _CryptDriftR_LinearGillespie(SEXP numSimSEXP, SEXP nuSEXP, SEXP a_matSEXP, SEXP x0SEXP, SEXP time_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numSim(numSimSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a_mat(a_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time_vec(time_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(LinearGillespie(numSim, nu, a_mat, x0, time_vec));
    return rcpp_result_gen;
END_RCPP
}
// Crypt_drift
arma::mat Crypt_drift(double alpha, double beta, int Ns, arma::vec time_intervals, bool persisting, int splitNum);
RcppExport SEXP _CryptDriftR_Crypt_drift(SEXP alphaSEXP, SEXP betaSEXP, SEXP NsSEXP, SEXP time_intervalsSEXP, SEXP persistingSEXP, SEXP splitNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time_intervals(time_intervalsSEXP);
    Rcpp::traits::input_parameter< bool >::type persisting(persistingSEXP);
    Rcpp::traits::input_parameter< int >::type splitNum(splitNumSEXP);
    rcpp_result_gen = Rcpp::wrap(Crypt_drift(alpha, beta, Ns, time_intervals, persisting, splitNum));
    return rcpp_result_gen;
END_RCPP
}
// Infer_Biased_mcmc
arma::mat Infer_Biased_mcmc(arma::umat x, arma::vec time_intervals, int max_iter, int Ns_WT, double lambda_WT, double tau_fix);
RcppExport SEXP _CryptDriftR_Infer_Biased_mcmc(SEXP xSEXP, SEXP time_intervalsSEXP, SEXP max_iterSEXP, SEXP Ns_WTSEXP, SEXP lambda_WTSEXP, SEXP tau_fixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time_intervals(time_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type Ns_WT(Ns_WTSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_WT(lambda_WTSEXP);
    Rcpp::traits::input_parameter< double >::type tau_fix(tau_fixSEXP);
    rcpp_result_gen = Rcpp::wrap(Infer_Biased_mcmc(x, time_intervals, max_iter, Ns_WT, lambda_WT, tau_fix));
    return rcpp_result_gen;
END_RCPP
}
// Infer_Biased_wLambda_mcmc
arma::mat Infer_Biased_wLambda_mcmc(arma::umat x, arma::vec time_intervals, int max_iter, int Ns_WT, double tau_fix);
RcppExport SEXP _CryptDriftR_Infer_Biased_wLambda_mcmc(SEXP xSEXP, SEXP time_intervalsSEXP, SEXP max_iterSEXP, SEXP Ns_WTSEXP, SEXP tau_fixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time_intervals(time_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type Ns_WT(Ns_WTSEXP);
    Rcpp::traits::input_parameter< double >::type tau_fix(tau_fixSEXP);
    rcpp_result_gen = Rcpp::wrap(Infer_Biased_wLambda_mcmc(x, time_intervals, max_iter, Ns_WT, tau_fix));
    return rcpp_result_gen;
END_RCPP
}
// Infer_Neutral_mcmc
arma::mat Infer_Neutral_mcmc(arma::umat x, arma::vec time_intervals, int max_iter);
RcppExport SEXP _CryptDriftR_Infer_Neutral_mcmc(SEXP xSEXP, SEXP time_intervalsSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type time_intervals(time_intervalsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Infer_Neutral_mcmc(x, time_intervals, max_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CryptDriftR_Biased_Cont_Labelling_Theory_c", (DL_FUNC) &_CryptDriftR_Biased_Cont_Labelling_Theory_c, 4},
    {"_CryptDriftR_LinearGillespie", (DL_FUNC) &_CryptDriftR_LinearGillespie, 5},
    {"_CryptDriftR_Crypt_drift", (DL_FUNC) &_CryptDriftR_Crypt_drift, 6},
    {"_CryptDriftR_Infer_Biased_mcmc", (DL_FUNC) &_CryptDriftR_Infer_Biased_mcmc, 6},
    {"_CryptDriftR_Infer_Biased_wLambda_mcmc", (DL_FUNC) &_CryptDriftR_Infer_Biased_wLambda_mcmc, 5},
    {"_CryptDriftR_Infer_Neutral_mcmc", (DL_FUNC) &_CryptDriftR_Infer_Neutral_mcmc, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CryptDriftR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}