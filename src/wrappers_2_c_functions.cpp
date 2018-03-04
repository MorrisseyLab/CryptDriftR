// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "biased_drift_functions.hpp"
#include "mcmcBiasedDrift.hpp"
#include "mcmcBiasedDrift_wLambda.hpp"
#include "mcmcNeutralDrift.hpp"
// using namespace arma;

// [[Rcpp::export]]
Rcpp::List Biased_Cont_Labelling_Theory_c(double kappa, double alpha, double beta, int Ns)
{
  double t;
  arma::rowvec theory_timeseries_row0;
  Rcpp::List ret;
  arma::vec    f_nm(Ns-1), g_nm(Ns-1), aux1(Ns-1), pow_gnm(Ns-1);
  double pi = arma::datum::pi;
  
  double slope_mono, intercept_mono;
  arma::vec partials(Ns-1);
  
  arma::vec              m  = arma::linspace<arma::vec>(1, Ns-1, Ns-1);
  double          cte = kappa*2./Ns;
  double  al_bet_sqrt   = sqrt(alpha*beta);
  double  beta_ov_alpha = beta/alpha;
  double  aux_al_bet    = alpha + beta - 2. * al_bet_sqrt;
  

  g_nm = 4*al_bet_sqrt*square(sin(0.5*pi*m/Ns)) + aux_al_bet;

  // .. For  1 <= n <= Ns-1
  for(int n = 1; n<Ns; n++)
  {
    f_nm = sin(pi*m/Ns)%sin(pi*m*n/Ns)/g_nm;
    partials(n-1) = arma::as_scalar(pow(beta_ov_alpha, (n-1)/2.)*cte*sum(f_nm));
  }

  // .. For n = Ns
  f_nm = sin(pi*m/Ns)%sin(pi*m*(Ns-1)/Ns);
  f_nm = f_nm/g_nm;
  
  slope_mono     =  beta*pow(beta_ov_alpha, (Ns-2)/2.)*cte*sum(f_nm);
  intercept_mono = -beta*pow(beta_ov_alpha, (Ns-2)/2.)*cte*sum(f_nm/g_nm);

  ret["Slope_mono"] = slope_mono;
  ret["Intercept_mono"] = intercept_mono;
  ret["Partials"] = partials;
  return ret;
}

// [[Rcpp::export]]
arma::umat LinearGillespie(int numSim, arma::umat nu, arma::mat a_mat, arma::uvec x0, arma::vec time_vec)
{
  int sim_i;
  arma::umat simulation_data_all=arma::zeros<arma::umat>(x0.n_elem, time_vec.n_elem);
  arma::umat simulation_data=arma::zeros<arma::umat>(x0.n_elem, time_vec.n_elem);
  for(sim_i = 0; sim_i<numSim; sim_i++)
  {
    LinearGillespie_indiv(nu, a_mat, x0, time_vec, simulation_data);
    simulation_data_all = simulation_data_all + simulation_data;
  }
  return(simulation_data_all);
}

// [[Rcpp::export]]
arma::mat Crypt_drift(double alpha, double beta, int Ns, arma::vec time_intervals, bool persisting, int splitNum = 0)
{
  arma::mat out_val;
  out_val = Crypt_drift_c(alpha, beta, Ns, time_intervals, persisting, splitNum);
  return(out_val);
}

// [[Rcpp::export]]
arma::mat Infer_Biased_mcmc(arma::umat x, arma::vec time_intervals, int max_iter, int Ns_WT, double lambda_WT, double tau_fix)
{
  arma::mat out_val;
  out_val = MH_pulse_chase_Biased(x, time_intervals, max_iter, Ns_WT, lambda_WT, tau_fix);
  return(out_val);
}

// [[Rcpp::export]]
arma::mat Infer_Biased_wLambda_mcmc(arma::umat x, arma::vec time_intervals, int max_iter, int Ns_WT, double tau_fix)
{
  arma::mat out_val;
  out_val = MH_pulse_chase_Biased_wLambda(x, time_intervals, max_iter, Ns_WT, tau_fix); 
  return(out_val);
}

// [[Rcpp::export]]
arma::mat Infer_Neutral_mcmc(arma::umat x, arma::vec time_intervals, int max_iter)
{
  arma::mat out_val;
  out_val = MH_pulse_chase_Neutral(x, time_intervals, max_iter);
  return(out_val);
}



