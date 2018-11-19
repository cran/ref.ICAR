#' Parameter Summaries for MCMC Analysis
#'
#' @description Takes a matrix of MCMC chains, iterations, and acceptance values
#' to return posterior summaries of the parameters, including posterior medians, intervals,
#' and acceptance rates.
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @importFrom stats median quantile
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda HPDinterval mcmc
#' @param MCMCchain Matrix of MCMC chains for the ICAR model parameters.
#' @param tauc.MCMC MCMC chains for the spatial dependence parameter.
#' @param sigma2.MCMC MCMC chains for the variance of the unstructured random effects.
#' @param beta.MCMC MCMC chains for the fixed effect regression coefficients.
#' @param phi.MCMC MCMC chains for the spatial random effects.
#' @param accept.phi Final acceptance number for spatial random effects.
#' @param accept.sigma2 Final acceptance number for variance of the unstructured random effects.
#' @param accept.tauc Final acceptance number for the spatial dependence parameter.
#' @param iters Number of MCMC iterations in \code{MCMCchain}.
#' @param burnin Number of MCMC iterations discarded as burn-in for \code{MCMCchain}.
#'
#' @return Parameter summaries
#'     \item{beta.median}{Posterior medians of the fixed effect regression coefficients.}
#'     \item{beta.hpd}{Highest Posterior Density intervals for the fixed effect regression coefficients.}
#'     \item{tauc.median}{Posterior median of the spatial dependence parameter.}
#'     \item{tauc.hpd}{Highest Posterior Density interval for the spatial dependence parameter.}
#'     \item{sigma2.median}{Posterior median of the unstructured random effects variance.}
#'     \item{sigma2.hpd}{Highest Posterior Density interval for the unstructured random effects variance.}
#'     \item{tauc.accept}{Final acceptance rate for the spatial dependence parameter.}
#'     \item{sigma2.accept}{Final acceptance rate for the unstructured random effects variance.}
#' @export
#' @examples
#' ## Refer to the vignette attached to the package.
ref.summary <- function(MCMCchain,tauc.MCMC,sigma2.MCMC,beta.MCMC,phi.MCMC,accept.phi,accept.sigma2,accept.tauc,iters=10000,burnin=5000) {

    beta.median <- apply(cbind(beta.MCMC[(burnin+1):iters,]),2,median)
    beta.obj <- mcmc(beta.MCMC,start=(burnin+1),end=iters)
    beta.hpd <- HPDinterval(beta.obj,prob=0.95)[,]
    tauc.median <- median(tauc.MCMC[(burnin+1):iters])
    tauc.obj <- mcmc(tauc.MCMC,start=(burnin+1),end=iters)
    tauc.hpd <- HPDinterval(tauc.obj,prob=0.95)[,]
    tauc.accept <- accept.tauc/iters
    sigma2.median <- median(sigma2.MCMC[(burnin+1):iters])
    sigma2.obj <- mcmc(sigma2.MCMC,start=(burnin+1),end=iters)
    sigma2.hpd <- HPDinterval(sigma2.obj,prob=0.95)[,]
    sigma2.accept <- accept.sigma2/iters

    return(list(beta.median=beta.median,beta.hpd=beta.hpd,tauc.median=tauc.median,tauc.hpd=tauc.hpd,
                sigma2.median=sigma2.median,sigma2.hpd=sigma2.hpd,
                tauc.accept=tauc.accept,sigma2.accept=sigma2.accept))
}
