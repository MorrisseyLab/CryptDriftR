% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DriftFunctions.R
\name{sim_contLabelling}
\alias{sim_contLabelling}
\title{Simulate the continuous labelling model using the gillespie algorithm. 
\code{sim_contLabelling}}
\usage{
sim_contLabelling(mu, lambda, Ns, time_points, numSim, Pr = 0.5, raw_sims = F)
}
\arguments{
\item{mu}{Stem cell mutation probability.}

\item{lambda}{Stem cell replacement rate.}

\item{Ns}{Number of functional stem cells per crypt.}

\item{time_points}{vector of time point values.}

\item{numSim}{Number of simulations to run.}

\item{Pr}{Bias in stem cell replacement Pr = 0.5 is neutral.}

\item{raw_sims}{If true returns raw simulation output. Otherwise just full and partial crypts (see description).}
}
\value{
If the raw_sims option was set to false (default), a matrix with two rows is returned the first being partial crypt frequency and the second monoclonal frequency. 
If raw_sims is set to true the matrix will have N + 1 rows, where the first row is the number of sims with no labelled crypts, the second 1 label crypt and so on.
}
\description{
The model can be found in the paper "Continuous clonal labeling reveals small numbers of functional stem cells in intestinal crypts and adenomas." 
Kozar, Morrissey et. al. 2013 Cell Stem Cell
}
\examples{
time_points = 100:200
x = sim_contLabelling(mu = 1.1e-4, lambda = 0.1, Ns = 5, time_points = 101:201, numSim = 1e7)


}
