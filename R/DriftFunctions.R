#' Simulate the continuous labelling model using the gillespie algorithm. 
#' \code{sim_contLabelling} 
#' 
#' The model can be found in the paper "Continuous clonal labeling reveals small numbers of functional stem cells in intestinal crypts and adenomas." 
#' Kozar, Morrissey et. al. 2013 Cell Stem Cell
#'
#' @param mu Stem cell mutation probability. 
#' @param Ns Number of functional stem cells per crypt.   
#' @param lambda Stem cell replacement rate.
#' @param time_points vector of time point values.
#' @param numSim Number of simulations to run.
#' @param Pr Bias in stem cell replacement Pr = 0.5 is neutral. 
#' @param raw_sims If true returns raw simulation output. Otherwise just full and partial crypts (see description).
#' 
#' @return If the raw_sims option was set to false (default), a matrix with two rows is returned the first being partial crypt frequency and the second monoclonal frequency. 
#' If raw_sims is set to true the matrix will have N + 1 rows, where the first row is the number of sims with no labelled crypts, the second 1 label crypt and so on. 
#' 
#' @examples
#' time_points = 100:200
#' x = sim_contLabelling(mu = 1.1e-4, lambda = 0.1, Ns = 5, time_points = 101:201, numSim = 1e7)
#' 
#'
#' @export
sim_contLabelling <- function(mu, lambda, Ns, time_points, numSim, Pr = 0.5, raw_sims = F)
{
  beta  <- 2*lambda*Pr
  alpha <- 2*lambda*(1-Pr)
  #   Number of chain reactions
  M      <- Ns
  # Initial state vector
  x0             <- c(numSim, rep(0,M)) 
  names(x0)      <- paste("x",seq(M+1),sep="") 
  nu             <- matrix(rep(0,(M*(M+1))),ncol=M)
  diag(nu)       <- -1
  diag(nu[2:M,]) <- +1
  nu[M+1,M]      <- +1
  nu             <- cbind(nu, -1*nu)
  nu             <- nu[,-1*c(ncol(nu))]
  a_mat          <- cbind(c(lambda*Ns*mu, rep(beta,M-1), rep(alpha,M-1)), c(seq(M),2:M))
  
  sims <- LinearGillespie(numSim=1, nu, a_mat, x0, time_points) 
  if(!raw_sims){
    return_val = rbind(colSums(sims[2:Ns,])/numSim, sims[Ns+1,]/numSim)
  }else{return_val = sims}
  return_val
}

#' Analytical solution to stem cell label drift model. 
#' \code{th_PulseChase}  Stochastic model used in "Defining stem cell dynamics in models of intestinal tumor initiation." 
#' Vermeulen, Morrissey, et. al. Science 2013. The model is a stochastic birth death model with absorbing boundaries at 0 and N. 
#'
#' @param lambda Stem cell replacement rate.
#' @param Ns Number of functional stem cells per crypt.   
#' @param tau Time delay until drift starts.   
#' @param time_points vector of time point values.
#' @param persisting Whether or not to return probability distribution of only the surviving clones or not.
#' @param Pr Bias in stem cell replacement Pr = 0.5 is neutral. 
#' @param splitNum To match to data choose the number of bins to distribute the model into (e.g 1/8s, splitNum = 8).  
#' @param alpha_beta two dim vector giving rate of replacement, alternative to using Pr. The relation is beta   = 2*lambda*Pr, alpha  = 2*lambda*(1-Pr) 
#' 
#' @return a matrix with the solution to the neutral or biased drift model. 
#' 
#' @examples
#' x = th_PulseChase(lambda = 0.1, Ns = 5, tau = 1, time_points = 1:50)
#' @export
th_PulseChase <- function(lambda, Ns, tau, time_points, persisting = T, Pr = 0.5, splitNum = 0, alpha_beta = NULL)
{
  if (time_points[1] < tau) stop("The first time point must be larger than tau!!!!")
  if(!is.null(alpha_beta))
  {
    alpha = alpha_beta[1]
    beta  = alpha_beta[2]
  }else{
    beta   = 2*lambda*Pr
    alpha  = 2*lambda*(1-Pr)    
  }
  Crypt_drift(alpha, beta, Ns, time_points - tau, persisting, splitNum)  
}


CA30Function_MonoClonal <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  intercept <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  aplha_lamb <- alpha*lambda
  for(i in 1:length(time.intervals))
  {
    t <- time.intervals[i]
    theory.timeseries[i] <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2)*(1-exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2)))
  }  
  aplha_lamb*time.intervals - theory.timeseries
}

CA30Function_PartialClones <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries   <- rep(0, length(time.intervals))
  theory.timeseries.n <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  n.all   <- 1:(Ns-1)
  for(n in n.all)
  {
    for(i in 1:length(time.intervals))
    {
      t                      <- time.intervals[i]
      theory.timeseries.n[i] <- 0.5*alpha*sum(sin(pi*m/Ns)*sin(pi*m*n/Ns)/(sin(0.5*pi*m/Ns))^2*exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2))
    }  
    theory.timeseries <- theory.timeseries+theory.timeseries.n   
  }
  # Return value
  alpha*Ns*(Ns-1)/2 - theory.timeseries
}

#' Analytical approximation to continuous labelling model. 
#'
#' \code{th_contLabelling} 
#' @param mu Stem cell replacement rate.
#' @param lambda Stem cell replacement rate.
#' @param Ns Number of functional stem cells per crypt.   
#' @param time_points vector of time point values.
#' 
#' @return a matrix with two rows, the frequency of partial clone cypts and the frequency of monoclonal crypts. 
#' 
#' @examples
#' x = th_contLabelling(mu = 1.1e-4, lambda = 0.1, Ns = 5, time_points = 101:201)
#'
#' @export
th_contLabelling <- function(mu, lambda, Ns, time_points)
{
  rbind(CA30Function_PartialClones(mu, lambda, Ns, time_points), 
        CA30Function_MonoClonal(mu, lambda, Ns, time_points))  
}
  
#' Analytical approximation to continuous labelling model, also allows bias. 
#' \code{th_contLabelling} 
#'
#'  The function returns only the slope of the monoclonal crypts and the partial clone levels. 
#'
#' @param mu Stem cell replacement rate.
#' @param lambda Stem cell replacement rate.
#' @param Ns Number of functional stem cells per crypt.   
#' @param Pr Bias in stem cell replacement Pr = 0.5 is neutral. 
#' @param alpha_beta two dim vector giving rate of replacement, alternative to using Pr. The relation is beta   = 2*lambda*Pr, alpha  = 2*lambda*(1-Pr) 
#' 
#' @return A list with the slope of the monoclonal accumulation, the intercept of the monoclonal accumulation and the level of the partials.
#' 
#' @examples
#' x = th_contLabelling_slope_and_partials(mu = 1.1e-4, lambda = 0.1, Ns = 5)
#'
#' @export
th_contLabelling_slope_and_partials <- function(mu, lambda, Ns, Pr = 0.5, alpha_beta = NULL)
{
  if(!is.null(alpha_beta))
  {
    alpha = alpha_beta[1]
    beta  = alpha_beta[2]
  }else{
    beta  = 2*lambda*Pr
    alpha = 2*lambda*(1-Pr)    
  }
  Biased_Cont_Labelling_Theory_c(mu, alpha, beta, Ns)
}


#' Run linear Gillespie simulations.
#' \code{LinearGillespie} 
#'
#' Use this function to run gillespie simulations. The function is used internally and is not strictly necessary to use directly. 
#' @param numSim Number of times to run the simulation. The output is the sum of all sims.
#' @param nu Matrix with a column per reaction and the same number of rows as states. Each column 
#' @param a_mat Two column matrix with the same number of rows as reactions. Each row contains first the rate of that reaction and second the   
#' @param x0 The initial state of the system. 
#' @param time_vec A vector of time points at which to report the simulation status.
#' 
#' Simulate a stochastic process with the Gillespie algorithm.    
#' 
#' @return A matrix with a column per data point and the same number of rows as states. 
#' 
#' @examples
#' 
#' # This is the code used to simulate continous labelling
#' mu     = 1e-4
#' lambda = 0.1
#' N      = 5
#' Pr     = 0.5
#' time_points = 100:150
#' numSim = 1e7
#' 
#' beta  <- 2*lambda*Pr
#' alpha <- 2*lambda*(1-Pr)
#' #   Number of chain reactions
#' M      <- N
#' # Initial state vector
#' x0             <- c(numSim, rep(0,M)) 
#' names(x0)      <- paste("x",seq(M+1),sep="") 
#' nu             <- matrix(rep(0,(M*(M+1))),ncol=M)
#' diag(nu)       <- -1
#' diag(nu[2:M,]) <- +1
#' nu[M+1,M]      <- +1
#' nu             <- cbind(nu, -1*nu)
#' nu             <- nu[,-1*c(ncol(nu))]
#' a_mat          <- cbind(c(lambda*N*mu, rep(beta,M-1), rep(alpha,M-1)), c(seq(M),2:M))
#' sims <- LinearGillespie(numSim=1, nu, a_mat, x0, time_points) 
#'
#' @export
LinearGillespie <- function(numSim, nu, a_mat, x0, time_vec) {
  LinearGillespie_c(numSim, nu, a_mat, x0, time_vec)
}

