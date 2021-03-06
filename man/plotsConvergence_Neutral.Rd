% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_NeutralDrift.R
\name{plotsConvergence_Neutral}
\alias{plotsConvergence_Neutral}
\title{MCMC convergence plots}
\usage{
plotsConvergence_Neutral(fit_vals)
}
\arguments{
\item{fit_vals}{list produced by fitNeutralDrift containing amongst other things the MCMC chains.}
}
\value{
Return ggplot object with convergence plots.
}
\description{
\code{plotsConvergence_Neutral} will plot convergence plots using the MCMC chains.
}
\examples{
time_points = c(4, 7, 10, 14, 21)
x = simulateNeutralDriftData(0.1, 5, 2, time_points, 8, 300)
fit_out = fitNeutralDrift(x, time_points)
plotsConvergence_Neutral(fit_out)
plotsNeutralDrift_Fit(fit_out)

}
