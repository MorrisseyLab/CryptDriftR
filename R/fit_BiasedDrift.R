#' Fit crypt biased stem cell drift model.
#'
#' \code{fitBiasedDrift} will run an mcmc algorithm to infer the parameters from the biased drift model. 
#' 
#' @param x Clone size counts. With a column per day and a row per size clone size.  
#' @param time_interval vector of time point values. Should have same number of elements as columns in x.
#' @param WT_params List with wild type neutral drift parameters (N, lambda and tau).
#' @param max_iter number of mcmc iterations to run.
#' @param n_par_chains number of parallel chains (and cpus used) to run. 
#' @param burn_in mcmc burn in 
#' @param thin mcmc thining 
#' 
#' @return Will retrun a list that can be used directly with the plotting functions. The list contains an mcmc matrix, the data and time vector used to fit, as well 
#' as the MAP estimate of the parameters.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' WT_params = list(lambda = 0.1, N  = 5, tau = 1)
#' x = simulateBiasedDriftData(WT_params$lambda, WT_params$N, WT_params$tau, Pr = 0.7, time_points, 8, 300)
#' 
#' fit_out = fitBiasedDrift(x, time_points, WT_params)
#' plotsConvergence_Biased(fit_out)
#' plotsBiasDrift_Fit(fit_out)
#'
#'@export
fitBiasedDrift = function(x, time_interval, WT_params, max_iter = 80000, n_par_chains = 4, burn_in = 5000, thin = 50)
{
  Ns.WT     = WT_params$N
  lambda.WT = WT_params$lambda
  tau.WT    = WT_params$tau
  x = as.matrix(x)
  registerDoMC(cores=min(detectCores(), n_par_chains))  
  all_runs = foreach(mcmc.run = 1:n_par_chains, .combine= rbind)%dopar%{
    mcmc_out = Infer_Biased_mcmc(x, time_interval, max_iter, Ns.WT, lambda.WT, tau.WT)
    indx_use = seq(from = burn_in, to = nrow(mcmc_out), by = thin)
    mcmc_sub = mcmc_out[indx_use, ]
    return(mcmc_sub)
  }  
  Pr_quant = getBiasParams(all_runs)
  list(mcmc = all_runs, x = x, time_interval = time_interval, WT_params = WT_params, Pr_quant = Pr_quant)
}

#' Simulate data from biased drift model for testing inference.
#'
#' \code{simulateBiasedDriftData} will use the analytical biased stem cell drift model solutions to similate count data. 
#' 
#' @param lambda Stem cell replacement rate.
#' @param Ns Number of functional stem cells per crypt.   
#' @param tau Time delay until drift starts. 
#' @param Pr Bias in stem cell replacement Pr = 0.5 is neutral. 
#' @param time_points vector of time point values.
#' @param size_measure Fractions in which the crypts are measured in (e.g. 8ths).
#' @param num_crypts Number of crypts measured per time point.
#' 
#' @return Will retrun simulated crypt clone count data.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' WT_params = list(lambda = 0.1, N  = 5, tau = 1)
#' x = simulateBiasedDriftData(WT_params$lambda, WT_params$N, WT_params$tau, Pr = 0.7, time_points, 8, 300)
#'
#'@export
simulateBiasedDriftData = function(lambda, Ns, tau, Pr, time_points, size_measure, num_crypts)
{
  Pn    <- th_PulseChase(lambda, Ns, tau, time_points, persisting = T, Pr = Pr, splitNum = size_measure)
  xn    <- matrix(0, size_measure, length(time_points))
  for(i in 1:length(time_points))
  {
    x.i          <- sample(1:size_measure, num_crypts, replace = T, prob = Pn[,i])
    t.i          <- table(x.i)
    indx.i       <- as.numeric(names(t.i))
    xn[indx.i,i] <- t.i    
  }
  xn
}

#' MCMC convergence plots
#'
#' \code{plotsConvergence_Biased} will plot convergence plots using the MCMC chains. 
#' 
#' @param fit_vals list produced by fitBiasedDrift containing amongst other things the MCMC chains.
#' 
#' @return Return ggplot object with convergence plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' WT_params = list(lambda = 0.1, N  = 5, tau = 1)
#' x = simulateBiasedDriftData(WT_params$lambda, WT_params$N, WT_params$tau, Pr = 0.7, time_points, 8, 300)
#' 
#' fit_out = fitBiasedDrift(x, time_points, WT_params)
#' plotsConvergence_Biased(fit_out)
#' plotsBiasDrift_Fit(fit_out)
#'
#'@export
plotsConvergence_Biased= function(fit_vals)
{
  
  sub.Pr        <- fit_vals$mcmc[,1]
  sub.unsc.post <- fit_vals$mcmc[,2]
  
  
  pp1 = qplot(y =        sub.Pr, x = 1:length(sub.Pr), geom = "line") + xlab("Iteration (thinned)") + ylab("Pr")
  pp2 = qplot(y = sub.unsc.post, x = 1:length(sub.Pr), geom = "line") + xlab("Iteration (thinned)") + ylab("Unscaled Posterior")
  
  pp3 = gg_acf(sub.Pr)
  pp4 = gg_acf(sub.unsc.post)
  
  pp5 = gg_qq(sub.Pr)
  pp6 = gg_qq(sub.unsc.post)
  
  
  grid.arrange(pp1, pp2, 
               pp3, pp4,
               pp5, pp6, nrow=3)
  
}


#' MCMC inference plots for biased drift model
#'
#' \code{plotPosterior_Pr} will plot posterior plots using the MCMC chains. 
#' 
#' @param fit_vals list produced by fitBiasedDrift containing amongst other things the MCMC chains.
#' 
#' @return Return ggplot object with posterior plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' WT_params = list(lambda = 0.1, N  = 5, tau = 1)
#' x = simulateBiasedDriftData(WT_params$lambda, WT_params$N, WT_params$tau, Pr = 0.7, time_points, 8, 300)
#' 
#' fit_out = fitBiasedDrift(x, time_points, WT_params)
#' plotsConvergence_Biased(fit_out)
#' plotPosterior_Pr(fit_out)
#'
#'@export
plotPosterior_Pr = function(fit_vals)
{
  all.runs.subs = data.frame(fit_vals$mcmc)
  all.mcmc   = data.frame(prob = all.runs.subs[,1])  
  vals.cols2 = c("#999999", "darkgreen", "steelblue" ,"#66CC00")
  pp =  ggplot(all.mcmc, aes(x = prob)) +
    geom_histogram(binwidth = 0.002, col = vals.cols2[3], fill = vals.cols2[3]) + 
    labs(x = "Mutant cell fitness", y = "Probability density") +
    geom_vline(xintercept = 0.5, lty =2, size = 1.5) + #ylim(0, 600)+
    ggtitle("Distribution for bias parameter") + xlim(0,1)
  pp
}


#' Plot data with biased drift model fit
#'
#' \code{plotsBiasDrift_Fit} will plot  data and fit of biased drift model. 
#' 
#' @param fit_vals list produced by fitBiasedDrift containing amongst other things the MCMC chains.
#' @param max_x Choose maximum value for the x axis (time)

#' @return Return ggplot object with posterior plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' WT_params = list(lambda = 0.1, N  = 5, tau = 1)
#' x = simulateBiasedDriftData(WT_params$lambda, WT_params$N, WT_params$tau, Pr = 0.7, time_points, 8, 300)
#' 
#' fit_out = fitBiasedDrift(x, WT_params, time_points)
#' plotsConvergence_Biased(fit_out)
#' plotsBiasDrift_Fit(fit_out)
#'
#'@export
plotsBiasDrift_Fit = function(fit_vals, max_x = 100)
{
  # Transform format for plotting
  uu2   = format_exp_data(data_x = as.data.frame(fit_vals$x), time_points = fit_vals$time_interval)
  WT_params = fit_vals$WT_params
  
  col_p = "steelblue"
  pp1 = ggplot(uu2, aes(y = value, x = Age)) + geom_point(size = 3, col = col_p) + 
    geom_errorbar(aes(ymin = low_lim, ymax = hi_lim), col = col_p) + 
    ggtitle("Biased drift fit") + 
    labs(x = "Time post labelling (days)", y = "Clone fraction") +
    xlim(0, max_x) + ylim(0,1) + theme_bw(24) +
    facet_grid(~ CloneSize)

  time_sim        = WT_params$tau:max_x
  
  sim_WT    = th_PulseChase(lambda = WT_params$lambda, Ns = WT_params$N, tau = WT_params$tau,
                                  time_points = time_sim, splitNum = nrow(fit_vals$x))
  sim_WT    = format_th_data(th_x = sim_WT, time_points = time_sim)
  

  sim_mut = th_PulseChase(lambda = WT_params$lambda, Ns = WT_params$N, tau = WT_params$tau,
                          Pr = fit_vals$Pr_quant[2],  time_points = time_sim, splitNum = nrow(fit_vals$x))
  sim_mut = format_th_data(th_x = sim_mut, time_points = time_sim)
  
  cols = c("WT"="orange","Mutant"=col_p, "Dummy = #111111")
  pp <- pp1 + scale_colour_manual(name="Type", values=cols)  +
    geom_line(data = sim_WT, aes(col = "WT"), lty = 2)   +
    geom_line(data = sim_mut, aes(col = "Mutant"))
  pp  
}



getBiasParams = function(all.runs.subs)
{
  sub.p        = all.runs.subs[,1]
  quantile(sub.p, probs = c(0.025, 0.5, 0.975))
}

