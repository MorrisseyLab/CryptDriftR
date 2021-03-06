#' Fit crypt neutral drift model.
#'
#' \code{fitNeutralDrift} will run an mcmc algorithm to infer the parameters from the neutral drift model. 
#' 
#' @param x Clone size counts. With a column per day and a row per size clone size.  
#' @param time_interval vector of time point values. Should have same number of elements as columns in x.
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
#' x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
#' fit_out = fitNeutralDrift(x, time_points)
#' plotsConvergence_Neutral(fit_out)
#' plotsNeutralDrift_Fit(fit_out)
#'
#'@export
fitNeutralDrift = function(x, time_interval, max_iter = 40000, n_par_chains = 2, burn_in = 5000, thin = 20)
{
  x = as.matrix(x)
  registerDoMC(cores=min(detectCores(), n_par_chains))  
  all_runs = foreach(mcmc_run = 1:n_par_chains, .combine= rbind)%dopar%{
    mcmc_out = Infer_Neutral_mcmc(x, time_interval, max_iter)
    indx_use = seq(from = burn_in, to = nrow(mcmc_out), by = thin)
    mcmc_sub = mcmc_out[indx_use, ]
    return(mcmc_sub)
    # cbind(mcmc_out$N, mcmc_out$lambda, mcmc_out$tau, mcmc_out$unsc.post)    
  }  
  MAP_est = getNeutralDirftParams(all_runs)
  list(mcmc = all_runs, x = x, time_interval = time_interval, MAP_est = MAP_est)
}


#' Simulate data from neutral drift model for testing inference.
#'
#' \code{simulateNeutralDriftData} will use the analytical neutral drift model solutions to similate count data. 
#' 
#' @param lambda Stem cell replacement rate.
#' @param Ns Number of functional stem cells per crypt.   
#' @param tau Time delay until drift starts.   
#' @param time_points vector of time point values.
#' @param size_measure Fractions in which the crypts are measured in (e.g. 8ths).
#' @param num_crypts Number of crypts measured per time point.
#' 
#' @return Will retrun simulated crypt clone count data.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
#'
#'@export
simulateNeutralDriftData <- function(lambda, Ns, tau, time_points, size_measure, num_crypts)
{
  Pn <- th_PulseChase(lambda, Ns, tau, time_points, persisting = T, Pr = 0.5, splitNum = size_measure)
  xn <- matrix(0, size_measure, length(time_points))
  for(i in 1:length(time_points))
  {
    x.i          <- sample(1:size_measure, num_crypts, replace = T, prob = Pn[,i])
    t.i          <- table(x.i)
    indx.i       <- as.numeric(names(t.i))
    xn[indx.i,i] <- t.i    
  }
  xn
}

gg_acf = function(chain)
{
  ac.p = acf(chain, plot=F)
  acd  = data.frame(lag=ac.p$lag, acf=ac.p$acf)
  pp   = ggplot(acd, aes(x=lag, y=acf)) + geom_bar(stat="identity") +
    geom_hline(yintercept=c(0.05, -0.05), linetype="dashed") 
  if(any(is.nan(ac.p$acf)))
  {
    pp = textGrob("Can't calculate acf")
  }
  pp
}

gg_qq = function(chain)
{
  len.1.chain = round(length(chain)/2)
  aux_df = as.data.frame(qqplot(chain[1:len.1.chain], chain[(len.1.chain+1):length(chain)], 
                                plot.it=FALSE))
  pp = ggplot(aux_df) + geom_point(aes(x=x, y=y)) + geom_abline(slope=1,colour="red") + xlab("First half of run") +
    ylab("Second half of run")
  pp
}

gg_qq_int = function(chain)
{
  len.1.chain = round(length(chain)/2)
  nn1 = as.vector(table(factor(chain[1:len.1.chain], levels = 2:30)))
  nn2 = as.vector(table(factor(chain[(len.1.chain+1):length(chain)], levels = 2:30)))
  
  aux_df = data.frame(x = nn1/sum(nn1), y = nn2/sum(nn2))
  pp = ggplot(aux_df) + geom_point(aes(x=x, y=y)) + geom_abline(slope=1,colour="red") + xlab("First half of run") +
    ylab("Second half of run")
  pp
}

#' MCMC convergence plots
#'
#' \code{plotsConvergence_Neutral} will plot convergence plots using the MCMC chains. 
#' 
#' @param fit_vals list produced by fitNeutralDrift containing amongst other things the MCMC chains.
#' 
#' @return Return ggplot object with convergence plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
#' fit_out = fitNeutralDrift(x, time_points)
#' plotsConvergence_Neutral(fit_out)
#' plotsNeutralDrift_Fit(fit_out)
#'
#'@export
plotsConvergence_Neutral = function(fit_vals)
{
  
  sub.Ns        <- fit_vals$mcmc[,1]
  sub.lambda    <- fit_vals$mcmc[,2]
  sub.tau       <- fit_vals$mcmc[,3]
  sub.unsc.post <- fit_vals$mcmc[,4]
  
  
  pp1 = qplot(y =        sub.Ns, x = 1:length(sub.Ns), geom = "line") + xlab("Iteration (thinned)") + ylab("Ns")
  pp2 = qplot(y =    sub.lambda, x = 1:length(sub.Ns), geom = "line") + xlab("Iteration (thinned)") + ylab("Lambda")
  pp3 = qplot(y =       sub.tau, x = 1:length(sub.Ns), geom = "line") + xlab("Iteration (thinned)") + ylab("Tau")
  pp4 = qplot(y = sub.unsc.post, x = 1:length(sub.Ns), geom = "line") + xlab("Iteration (thinned)") + ylab("Unscaled Posterior")
  
  pp5 = gg_acf(sub.Ns)
  pp6 = gg_acf(sub.lambda)
  pp7 = gg_acf(sub.tau)
  pp8 = gg_acf(sub.unsc.post)
 
  pp9  = gg_qq_int(sub.Ns) #+ position_jitter()
  pp10 = gg_qq(sub.lambda)
  pp11 = gg_qq(sub.tau)
  pp12 = gg_qq(sub.unsc.post)
  
  
  arrangeGrob(pp1, pp2, pp3, pp4,
               pp5, pp6, pp7, pp8,
               pp9, pp10, pp11, pp12, nrow=3)
  
}

#' MCMC inference plots for neutral drift model
#'
#' \code{plotPosterior_Neutral} will plot posterior plots using the MCMC chains. 
#' 
#' @param fit_vals list produced by fitNeutralDrift containing amongst other things the MCMC chains.
#' 
#' @return Return ggplot object with posterior plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
#' fit_out = fitNeutralDrift(x, time_points)
#' plotsConvergence_Neutral(fit_out)
#' plotsNeutralDrift_Fit(fit_out)
#'
#'@export
plotPosterior_Neutral = function(fit_vals)
{
  all.runs.subs = data.frame(fit_vals$mcmc)
  colnames(all.runs.subs) = c("StemCells", "lambda", "tau", "UnscPost")
  ## Heat map
  n.vals       = 2:16
  breaks.hist  = seq(from = 0, to = 1.2*max(all.runs.subs[,2]), length.out=500)

  dist_2d = makePost(all.runs.subs, n.vals, breaks.hist)
  pp1     = ggplot(data =  dist_2d, aes(y = SC, x = RepRate)) + 
    geom_tile(aes(fill = Post)) +
    scale_fill_gradient(low = "white", high = "steelblue") + 
    xlab("Stem cell replacement rate (per day per stem cell)") + 
    ylab("Number of Stem Cells") +  xlim(0,0.4) +theme_bw()
#   print(p)

  all.runs.subs$StemCells =  factor(all.runs.subs$StemCells) 
  ## Marginals
  pp2 = ggplot(all.runs.subs, aes(x = StemCells)) + 
    stat_count(aes(y=..count../sum(..count..))) +  scale_x_discrete(limits = 1:8 )+ theme_bw(base_size = 12)
  pp3 = ggplot(all.runs.subs, aes(x = lambda)) + 
    geom_histogram(aes(y=..density..), binwidth=0.0015)+ theme_bw(base_size = 12)
  pp4 = ggplot(all.runs.subs, aes(x = tau)) + 
    geom_histogram(aes(y=..density..), binwidth=0.05)+ theme_bw(base_size = 12)
  arrangeGrob(pp1, pp2, pp3, pp4, nrow=2, top = textGrob("Neutral drift fit parameters", gp=gpar(fontsize=24), just="top"))
}


post.2D = function(sub.Ns, sub.lambda, bins = 300, n.vals = 2:30, breaks.hist = NULL)
{
  #   n.vals   <- 2:30
  max.lambda  <- max(sub.lambda)
  if(is.null(breaks.hist)) breaks.hist <-  seq(0.01, max.lambda, length.out = bins + 1)
  post.mat    <- matrix(0, length(n.vals), length(breaks.hist)-1)
  for (i in n.vals)
  {
    indx.n          <- sub.Ns == i
    hist.i          <- hist(sub.lambda[indx.n], breaks = breaks.hist, plot = F)
    post.mat[i-1, ] <- hist.i$counts
  }
  post.mat           = post.mat/length(sub.Ns)
  rownames(post.mat) = n.vals #paste("N", n.vals)
  colnames(post.mat) = round(hist.i$mids, 4)
  post.mat
}

makePost = function(all.runs.subs, n.vals, breaks.hist)
{
  sub.Ns        = all.runs.subs[,1]
  sub.lambda    = all.runs.subs[,2]

  ## Get post 2D mat
  post.mat = post.2D(sub.Ns, sub.lambda, bins = NA, n.vals, breaks.hist)
  post.mat = data.frame(SC = rownames(post.mat), post.mat)
  post.mat = post.mat %>% 
             gather(RepRate, Post, -SC) %>%
             mutate(SC = factor(SC, levels = n.vals)) %>%
             mutate(RepRate = str_replace(RepRate, "X", ""))%>%
             mutate(RepRate = as.numeric(str_replace(RepRate, fixed("e."), "e-")))
    
  post.mat$SC <- factor(post.mat$SC, levels = n.vals)
  
  ## Scale by max 
  #post.mat$Post <- post.mat$Post/max(post.mat$Post)  
  post.mat
}


#' Plot data with neutral drift model fit
#'
#' \code{plotsNeutralDrift_Fit} will plot  data and fit of neutral drift model. 
#' 
#' @param fit_vals list produced by fitNeutralDrift containing amongst other things the MCMC chains.
#' @param max_x Choose maximum value for the x axis (time)
#' 
#' @return Return ggplot object with posterior plots.
#' 
#' @examples
#' time_points = c(4, 7, 10, 14, 21)
#' x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
#' fit_out = fitNeutralDrift(x, time_points)
#' plotsConvergence_Neutral(fit_out)
#' plotsNeutralDrift_Fit(fit_out)
#'
#'@export
plotsNeutralDrift_Fit = function(fit_vals, max_x = 100)
{

  WT_params = fit_vals$MAP_est
  # Transform format for plotting
  
  uu2   = format_exp_data(data_x = as.data.frame(fit_vals$x), time_points = fit_vals$time_interval)
  col_p = "steelblue"
  pp1 = ggplot(uu2, aes(y = value, x = Age)) + geom_point(size = 3, col = col_p) + 
    geom_errorbar(aes(ymin = low_lim, ymax = hi_lim), col = col_p) + 
    ggtitle("Neutral drift fit") + 
    labs(x = "Time post labelling (days)", y = "Clone fraction") +
    xlim(0, max_x) + ylim(0,1) + theme_bw(24) +
    facet_grid(~ CloneSize)
  
  time_sim        = WT_params$tau:max_x
  sim_WT          = th_PulseChase(lambda = WT_params$lambda, Ns = WT_params$N, tau = WT_params$tau,
                                  time_points = time_sim, splitNum = nrow(fit_vals$x))
  sim_WT    = format_th_data(th_x = sim_WT, time_points = time_sim)
  pp = pp1 + geom_line(data = sim_WT) #, col = col_p)
  pp  
}

#' MAP estimate of parameters
getNeutralDirftParams = function(fit_out)
{
  params = fit_out[which.max(fit_out[,4]),1:3]
  list(lambda = params[2], N = params[1], tau = params[3])
}

