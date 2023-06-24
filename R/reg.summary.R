#' Regional Summaries for Areal Data Modeled by ICAR Reference Prior Model
#' @importFrom stats median quantile
#' @importFrom coda HPDinterval mcmc
#' @description This function takes data and sampled MCMC chains for an areal data set
#' and gives fitted posterior values and summaries by region using the model by \insertCite{keefe2018}{ref.ICAR}.
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @param MCMCchain Matrix of MCMC chains, using the sampling from \insertCite{keefe2018}{ref.ICAR}.
#' @param X Matrix of covariates.
#' @param Y Vector of responses.
#' @param burnin Number of MCMC iterations discarded as burn-in in \code{MCMCchain}.
#'
#' @return A list of the fitted distributions by region, and medians and credible intervals by region.
#'     \item{fit.dist}{Matrix of fitted posterior values for each region in the data.}
#'     \item{reg.medians}{Vector of posterior medians for fitted response by region.}
#'     \item{reg.cred}{Data frame of credbile intervals by region.}
#' @export
#' @examples
#' ## Refer to the vignette attached to the package.
reg.summary <- function(MCMCchain,X,Y,burnin) {
  num.reg <- length(Y)

  # specify fixed effect regression coefficients
  betavec <- MCMCchain[((burnin+1):nrow(MCMCchain)),3:(3+ncol(X)-1)]

  # specify model values for spatial effects
  phi.composite <- MCMCchain[((burnin+1):nrow(MCMCchain)),(3+ncol(X)):ncol(MCMCchain)]

  # fitted posteriors by region
  fit.dist <- matrix(0,nrow=nrow(MCMCchain)-burnin, ncol=num.reg)
  for(k in 1:(nrow(MCMCchain)-burnin)) {
    fit.dist[k,] <- X%*%betavec[k,] + as.matrix(phi.composite[k,])
  }

  # regional medians and intervals
  reg.medians <- apply(fit.dist,2,median)

  reg.cred <- data.frame(t(apply(fit.dist, 2, FUN = quantile, prob = c(0.025, 0.975))))

  reg.obj <- mcmc(fit.dist,start=(burnin+1),end=nrow(fit.dist))
  reg.hpd <- HPDinterval(reg.obj,prob=0.95)[,]

  return(list(fit.dist=fit.dist,reg.medians=reg.medians,reg.hpd=reg.hpd))
}
