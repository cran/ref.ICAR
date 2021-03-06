% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reg.summary.R
\name{reg.summary}
\alias{reg.summary}
\title{Regional Summaries for Areal Data Modeled by ICAR Reference Prior Model}
\usage{
reg.summary(MCMCchain, X, Y, burnin)
}
\arguments{
\item{MCMCchain}{Matrix of MCMC chains, using the sampling from \insertCite{keefe2018}{ref.ICAR}.}

\item{X}{Matrix of covariates.}

\item{Y}{Vector of responses.}

\item{burnin}{Number of MCMC iterations discarded as burn-in in \code{MCMCchain}.}
}
\value{
A list of the fitted distributions by region, and medians and credible intervals by region.
    \item{fit.dist}{Matrix of fitted posterior values for each region in the data.}
    \item{reg.medians}{Vector of posterior medians for fitted response by region.}
    \item{reg.cred}{Data frame of credbile intervals by region.}
}
\description{
This function takes data and sampled MCMC chains for an areal data set
and gives fitted posterior values and summaries by region using the model by \insertCite{keefe2018}{ref.ICAR}.
}
\examples{
## Refer to the vignette attached to the package.
}
\author{
Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
}
