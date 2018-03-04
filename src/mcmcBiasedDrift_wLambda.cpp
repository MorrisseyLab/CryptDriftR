// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "common_mcmc_funcs.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


mat MH_pulse_chase_Biased_wLambda(umat x, vec time_intervals, int max_iter, int Ns_WT, double tau_fix)
{
  double p_new, lambda_new, log_unsc_post, log_unsc_post_new;
  double sd_p, sd_lambda;
  mat Sigma_mat, chol_sigma;
  vec param_lik(2), param2_lik(3), sd_vals, new_p_lamb(2), p_lamb(2);
  int prev_iter;
  int p_indx      = 0;
  int lambda_indx = 1;
  int LogLik_indx = 2;
  mat mcmc_out(max_iter, 3);
  sd_vals << 0.01 << 0.25 << endr;
  
  // For the joint update 
  sd_p        = 0.1;
  sd_lambda   = 0.02;
  
  Sigma_mat << pow(sd_p,2) << -0.9*sd_lambda*sd_p << endr << -0.9*sd_lambda*sd_p << pow(sd_lambda, 2)<<endr;
  chol_sigma = chol(Sigma_mat).t();

  mcmc_out(0,p_indx)      = R::rbeta(0.5,0.5);
  mcmc_out(0,lambda_indx) = R::rbeta(0.5,0.5);
  
  log_unsc_post = unscaled_log_posterior_bias(x, time_intervals, mcmc_out(0,lambda_indx), Ns_WT, mcmc_out(0,p_indx), tau_fix);    
  mcmc_out(0,LogLik_indx) = log_unsc_post;
  for(int iter_i = 1; iter_i < max_iter;  iter_i++)
  {
    prev_iter = iter_i-1;
    // Joint update for p and lambda =================================
    p_lamb(0) = mcmc_out(prev_iter, p_indx);
    p_lamb(1) = mcmc_out(prev_iter, lambda_indx);
    new_p_lamb = joint_lamb_p_RW(p_lamb, chol_sigma); 
    // Unpack values
    p_new = new_p_lamb(0); lambda_new = new_p_lamb(1); 
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, lambda_new, Ns_WT, p_new, tau_fix);
    // Rcpp::Rcout<<log_unsc_post_new;
    
    param2_lik         = MH_double_update_step(p_lamb, new_p_lamb, log_unsc_post, log_unsc_post_new); 
    mcmc_out(iter_i, p_indx)      = param2_lik(0);
    mcmc_out(iter_i, lambda_indx) = param2_lik(1);
    log_unsc_post = param2_lik(2);
    // Just update lambda 
    lambda_new        = RW_propose(mcmc_out(iter_i, lambda_indx), sd_vals);
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, lambda_new, Ns_WT, mcmc_out(iter_i, p_indx), tau_fix);
    param_lik         = MH_update_step(mcmc_out(iter_i, lambda_indx), lambda_new, log_unsc_post, log_unsc_post_new);
    // Unpack values
    mcmc_out(iter_i, lambda_indx) = param_lik(0); log_unsc_post = param_lik(1);
    mcmc_out(iter_i, LogLik_indx) = log_unsc_post;
  }  
  return(mcmc_out);
}


