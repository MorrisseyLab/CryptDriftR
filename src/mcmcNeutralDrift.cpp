// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "common_mcmc_funcs.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

mat MH_pulse_chase_Neutral(umat x, vec time_intervals, int max_iter)
{
  double              lambda_new, tau_new, log_unsc_post;
  int                                              N_new;
  double    log_unsc_post_new, sd_RW, log_ratio_proposal;

  vec param_lik(2), new_N_lamb(2), N_lamb(2), new_vals(3), param2_lik(3);
  vec sd_vals_tau(2), sd_vals_lamb(2);
  int prev_iter;
  int N_indx      = 0;
  int lambda_indx = 1;
  int tau_indx    = 2;
  int LogLik_indx = 3;
  mat mcmc_out(max_iter, 4);
  sd_vals_tau  <<    1 <<    3 << endr;
  sd_vals_lamb << 0.01 << 0.25 << endr;
  sd_RW = 8e-5;
  // Initialise mcmc
  mcmc_out(0,      N_indx) = std::floor(std::abs(R::rnorm(7,3))) + 1;
  mcmc_out(0, lambda_indx) =                       R::rbeta(0.5,0.5);
  mcmc_out(0,    tau_indx) =                                      2.;
  
  log_unsc_post = unscaled_log_posterior_bias(x, time_intervals, mcmc_out(0,lambda_indx), mcmc_out(0,N_indx), 0.5, mcmc_out(0,tau_indx));    
  mcmc_out(0, LogLik_indx) = log_unsc_post;
  
  for(int iter_i = 1; iter_i < max_iter;  iter_i++)
  {
    prev_iter = iter_i-1;
    // Joint update for N and lambda =================================
    // Pack into vector
    N_lamb(0) = mcmc_out(prev_iter, N_indx); N_lamb(1) = mcmc_out(prev_iter,lambda_indx);
    // Propose new values and logratio of the proposal distributions
    new_vals  = joint_lamb_N_RW(mcmc_out(prev_iter, N_indx), mcmc_out(prev_iter,lambda_indx), sd_RW);
    // Unpack
    N_new     = new_vals(0); lambda_new = new_vals(1); log_ratio_proposal = new_vals(2);
    //Pack into vector 
    new_N_lamb(0) = N_new; new_N_lamb(1) = lambda_new;
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, lambda_new, N_new, 0.5, mcmc_out(prev_iter,tau_indx));
    // Update variables
    param2_lik = MH_double_update_step_non_symm(N_lamb, new_N_lamb, log_unsc_post, log_unsc_post_new, log_ratio_proposal); 
    // Unpack
    mcmc_out(iter_i, N_indx)      = param2_lik(0);
    mcmc_out(iter_i, lambda_indx) = param2_lik(1);
    log_unsc_post                 = param2_lik(2);
    // Just update lambda ==============================================
    lambda_new        = RW_propose(mcmc_out(iter_i, lambda_indx), sd_vals_lamb);
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, lambda_new, mcmc_out(iter_i, N_indx), 0.5, mcmc_out(prev_iter,tau_indx));
    param_lik         = MH_update_step(mcmc_out(iter_i, lambda_indx), lambda_new, log_unsc_post, log_unsc_post_new);
    // Unpack values
    mcmc_out(iter_i, lambda_indx) = param_lik(0); log_unsc_post = param_lik(1);
    // Just update tau ==============================================
    tau_new           = RW_propose(mcmc_out(prev_iter, tau_indx), sd_vals_tau);
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, mcmc_out(iter_i, lambda_indx), mcmc_out(iter_i, N_indx), 0.5, tau_new);
    param_lik         = MH_update_step(mcmc_out(prev_iter, tau_indx), tau_new, log_unsc_post, log_unsc_post_new);
    mcmc_out(iter_i, tau_indx) = param_lik(0); log_unsc_post = param_lik(1);
    // Store unsacles post
    mcmc_out(iter_i, LogLik_indx) = log_unsc_post;
  }  
  return(mcmc_out);
}


