#' Trace Plots for Parameters in ICAR Model
#' @importFrom graphics plot
#' @description This function creates trace plots for the parameters in
#' the ICAR reference prior model \insertCite{keefe2018}{ref.ICAR}.
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @param MCMCchain Matrix of MCMC chains for the model parameters.
#' @param X Matrix of covariates.
#' @param burnin Number of MCMC iterations from \code{MCMCchain} discarded as burn-in.
#' @param num.reg Number of regions in the areal data set.
#'
#' @return Trace plots for the fixed effect regression coefficients, the precision parameter,
#' and the unstructured random effects variance.
#' @export
#' @examples
#' ## Refer to the vignette attached to the package.
ref.plot <- function(MCMCchain, X, burnin, num.reg) {
    plot(MCMCchain[(burnin+1):nrow(MCMCchain),1], type='l', ylab="tauc")
    plot(MCMCchain[(burnin+1):nrow(MCMCchain),2], type='l', ylab="sigma2")
    #check for an intercept to correctly label plots
    if (mean(as.matrix(X[ ,1]))==1) {
        for (i in 3:(dim(MCMCchain)[2]-num.reg)) {
            plot(MCMCchain[(burnin+1):nrow(MCMCchain),i], type='l', ylab=paste("beta ",i-3))}
    }else {for (i in 3:(dim(MCMCchain)[2]-num.reg)) {
            plot(MCMCchain[(burnin+1):nrow(MCMCchain),i], type='l', ylab=paste("beta ",i-2))}
    }
}
