// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "common_mcmc_funcs.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

mat MH_pulse_chase_Biased(umat x, vec time_intervals, int max_iter, int Ns_WT, double lambda_WT, double tau_fix)
{
  double p_new, log_unsc_post, log_unsc_post_new, h_ratio, rnum_i;
  vec parma_lik, sd_vals;
  int prev_iter;
  int p_indx =0;
  int LogLik_indx =1;
  mat mcmc_out(max_iter, 2);
  sd_vals << 0.01 << 0.25 << endr;
  
  
  mcmc_out(0,p_indx) = R::rbeta(0.5,0.5);
  log_unsc_post = unscaled_log_posterior_bias(x, time_intervals, lambda_WT, Ns_WT, mcmc_out(0,p_indx), tau_fix);    
  mcmc_out(0,1) = log_unsc_post;
  for(int iter_i = 1; iter_i < max_iter;  iter_i++)
  {
    prev_iter = iter_i-1;
    // Update p =================================
    p_new             = RW_propose(mcmc_out(prev_iter, p_indx), sd_vals);
    log_unsc_post_new = unscaled_log_posterior_bias(x, time_intervals, lambda_WT, Ns_WT, p_new, tau_fix);
    parma_lik         = MH_update_step(mcmc_out(prev_iter, p_indx), p_new, log_unsc_post, log_unsc_post_new);
    mcmc_out(iter_i, p_indx) = parma_lik(0); log_unsc_post = parma_lik(1); // Unpack
    mcmc_out(iter_i, LogLik_indx) = log_unsc_post;
  }  
  return(mcmc_out);
}

