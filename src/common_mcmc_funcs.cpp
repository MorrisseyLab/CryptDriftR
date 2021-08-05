// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <cmath>
#include "biased_drift_functions.hpp"
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// Use this function for all pulse chase inference. 
// Although in some cases some of the parameters are not inferred, the priors will just multiply by a constant term.
double unscaled_log_posterior_bias(umat x, vec time_intervals, double lambda, int Ns, double p, double tau)
{
  int N_limit = 30; // Set limit on the value of N
  double unsc_log_post, size_measure;
  double shape_val, rate_val, aux_tau, log_prior, log_lik, aux_lamb;
  double beta   = 2*lambda*p;
  double alpha  = 2*lambda*(1-p);
  mat log_all_preds, probs_drift;
  
  if(lambda<0 | Ns<2 | Ns > N_limit | tau > time_intervals(0) | tau < 0 | p > 1 | p < 0)
  {
    unsc_log_post = -999999999999;
  }else{
    size_measure  = x.n_rows;
    probs_drift = Crypt_drift_c(alpha, beta, Ns, time_intervals-tau, true, size_measure);
    // If prob is zero make 1 so as to avoid NaN/Inf when taking log
    // probs_drift.elem(find(probs_drift < 1e-15)).fill(1);
    log_all_preds = log(probs_drift);
    log_all_preds.elem(find_nonfinite(log_all_preds)).fill(0);
    log_lik       = accu(log_all_preds%x);
    shape_val     = 0.0001; 
    rate_val      = 0.0001; 
    // aux_tau       = R::dgamma(   tau, shape_val, 1./rate_val, true);
    // aux_lamb      = R::dgamma(lambda, shape_val, 1./rate_val, true);
    aux_tau       = R::dnorm4(    tau,   0,   5, true); // Half normal
    aux_lamb      = R::dnorm4( lambda,   0,   5, true); // Half normal
    log_prior     = R::dbeta(       p, 0.5, 0.5, true) + aux_tau + aux_lamb;
    // N has a constant prior
    unsc_log_post = log_lik + log_prior;
  }
  return(unsc_log_post);
}


vec MH_update_step(double old_param, double new_param, double old_log_lik, double new_log_lik)
{
  vec param_and_loglik_out(2);
  double h_ratio, rnum;
  h_ratio           = std::min(0., new_log_lik - old_log_lik);
  rnum            = std::log(R::runif(0, 1));
  if(rnum < h_ratio)
  {
    param_and_loglik_out(0) = new_param;
    param_and_loglik_out(1) = new_log_lik;
  }
  else{
    param_and_loglik_out(0) = old_param;
    param_and_loglik_out(1) = old_log_lik;
  }
  return(param_and_loglik_out);
}

double RW_propose(double param_old, vec sd_vals)
{
  double param_new;
  vec prob_choose_vals, sd_i;
  prob_choose_vals << 0.8 << 0.2 << endr;
  sd_i      = Rcpp::RcppArmadillo::sample_main(sd_vals, 1, false, prob_choose_vals);
  param_new = R::rnorm(param_old, sd_i(0));
  return(param_new);
}

vec MH_double_update_step(vec old_param, vec new_param, double old_log_lik, double new_log_lik)
{
  vec param_and_loglik_out(3);
  double h_ratio, rnum;
  h_ratio  = std::min(0., new_log_lik - old_log_lik);
  rnum     = std::log(R::runif(0, 1));
  if(rnum < h_ratio)
  {
    param_and_loglik_out.subvec(0,1) = new_param;
    param_and_loglik_out(2) = new_log_lik;
  }
  else{
    param_and_loglik_out.subvec(0,1) = old_param;
    param_and_loglik_out(2) = old_log_lik;
  }
  return(param_and_loglik_out);
}

vec joint_lamb_p_RW(vec curr_p_lambda, mat Sigma_mat_chol)
{
  vec  outvec, rvec;
  rvec = Rcpp::rnorm(2, 0., 1.);
  outvec = curr_p_lambda + Sigma_mat_chol*rvec;
  return(outvec);
}

vec joint_lamb_N_RW(int curr_N, double curr_lambda, double sd_RW)
{
  double  new_lambda, curr_s, log_ratio_proposal;
  double     q_old_g_new, q_new_g_old, new_s;
  int                           new_N;
  vec                     out_vals(3);
  uvec                  choose_one(2);
  choose_one(0) = curr_N - 1;
  choose_one(1) = curr_N + 1;
  // Propose new N
  new_N       = Rcpp::RcppArmadillo::sample(choose_one, 1, false)(0);
  curr_s      = curr_lambda/pow(curr_N, 2);
  new_lambda  = R::rnorm(curr_s*pow(new_N,2), sd_RW*pow(new_N,2));
  new_s       = new_lambda/pow(new_N, 2);
  // Rcpp::Rcout << curr_N << " " << curr_lambda << " " << curr_s << "  " << new_N << "  " << new_lambda << " " << new_s << std::endl;
  q_old_g_new = R::dnorm4(curr_lambda,  new_s*pow(curr_N, 2),  sd_RW*pow(curr_N, 2), false);
  q_new_g_old = R::dnorm4( new_lambda,  curr_s*pow( new_N, 2), sd_RW*pow( new_N, 2), false);
  log_ratio_proposal = log(q_old_g_new/q_new_g_old);
  out_vals(0) = new_N; out_vals(1) = new_lambda; out_vals(2) = log_ratio_proposal;
  return(out_vals);
}

vec MH_double_update_step_non_symm(vec old_param, vec new_param, double old_log_lik, double new_log_lik, double log_ratio_proposal)
{
  vec param_and_loglik_out(3);
  double h_ratio, rnum;
  h_ratio  = std::min(0., new_log_lik - old_log_lik + log_ratio_proposal);
  rnum     = std::log(R::runif(0, 1));
  if(rnum < h_ratio)
  {
    param_and_loglik_out.subvec(0,1) = new_param;
    param_and_loglik_out(2) = new_log_lik;
  }
  else{
    param_and_loglik_out.subvec(0,1) = old_param;
    param_and_loglik_out(2) = old_log_lik;
  }
  return(param_and_loglik_out);
}


