double unscaled_log_posterior_bias(arma::umat x, arma::vec time_intervals, double lambda, int Ns, double p, double tau);
double RW_propose(double param_old, arma::vec sd_vals);
arma::vec MH_update_step(double old_param, double new_param, double old_log_lik, double new_log_lik);
arma::vec joint_lamb_p_RW(arma::vec curr_p_lambda, arma::mat Sigma_mat_chol);
arma::vec MH_double_update_step(arma::vec old_param, arma::vec new_param, double old_log_lik, double new_log_lik);
arma::vec joint_lamb_N_RW(int curr_N, double curr_lambda, double sd_RW);
arma::vec MH_double_update_step_non_symm(arma::vec old_param, arma::vec new_param, double old_log_lik, double new_log_lik, double log_ratio_proposal);
  

