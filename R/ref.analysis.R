#' MCMC Analysis and Summaries for Reference Prior on an Intrinsic Autoregressive Model for Areal Data
#'
#' @description Performs analysis on a geographical areal data set using the objective prior
#' for intrinsic conditional autoregressive (ICAR) random effects
#' \insertCite{keefe2018}{ref.ICAR}.  It takes a shapefile, data, and region names to
#' construct a neighborhood matrix and perform Markov chain Monte Carlo sampling on the
#' unstructured and spatial random effects.  Finally, the function obtains regional estimates and
#' performs posterior inference on the model parameters.
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#' @importFrom MCMCglmm posterior.mode
#' @importFrom coda HPDinterval mcmc
#' @importFrom sf st_read
#' @importFrom spdep poly2nb
#' @importFrom mvtnorm rmvnorm
#' @param shape.file A shapefile corresponding to the regions for analysis.
#' @param X A matrix of covariates, which should include a column of 1's for models with a non-zero intercept
#' @param y A vector of responses.
#' @param x.reg.names A vector specifying the order of region names contained in \code{X}.
#' @param y.reg.names A vector specifying the order of region names contained in \code{y}.
#' @param shp.reg.names A vector specifying the order of region names contained in the shapefile, if there is not a NAME column in the file.
#' @param iters Number of MCMC iterations to perform.  Defaults to 10,000.
#' @param burnin Number of MCMC iterations to discard as burn-in.  Defaults to 5,000.
#' @param verbose If FALSE, MCMC progress is not printed.
#' @param tauc.start Starting MCMC value for the spatial dependence parameter.
#' @param beta.start Starting MCMC value for the fixed effect regression coefficients.
#' @param sigma2.start Starting MCMC value for the variance of the unstructured random effects.
#' @param step.tauc Step size for the spatial dependence parameter.
#' @param step.sigma2 Step size for the variance of the unstructured random effects.
#'
#' @return A list containing \code{H}, MCMC chains, parameter summaries, fitted regional values,
#' and regional summaries.
#'     \item{H}{The neighborhood matrix.}
#'     \item{MCMC}{Matrix of MCMC chains for all model parameters.}
#'     \item{beta.median}{Posterior medians of the fixed effect regression coefficients.}
#'     \item{beta.hpd}{Highest Posterior Density intervals for the fixed effect regression coefficients.}
#'     \item{tauc.median}{Posterior median of the spatial dependence parameter.}
#'     \item{tauc.hpd}{Highest Posterior Density interval for the spatial dependence parameter.}
#'     \item{sigma2.median}{Posterior median of the unstructured random effects variance.}
#'     \item{sigma2.hpd}{Highest Posterior Density interval for the unstructured random effects variance.}
#'     \item{tauc.accept}{Final acceptance rate for the spatial dependence parameter.}
#'     \item{sigma2.accept}{Final acceptance rate for the unstructured random effects variance.}
#'     \item{fit.dist}{Matrix of fitted posterior values for each region in the data.}
#'     \item{reg.medians}{Vector of posterior medians for fitted response by region.}
#'     \item{reg.hpd}{Data frame of Highest Posterior Density intervals by region.}
#'
#' @export
#' @examples
#' ## Refer to the vignette attached to the package.
#'
ref.analysis <- function(shape.file,X,y,x.reg.names,y.reg.names,shp.reg.names=NULL,iters=10000,
                         burnin=5000,verbose=TRUE,tauc.start=1,beta.start=1,sigma2.start=1,step.tauc=0.5,step.sigma2=0.5) {

  ## Check for missing data
  if(mean(is.na(y))!=0) {stop("The current version of the package does not support missing data.")}
  num.reg <- nrow(as.matrix(y))

  ## Get neighborhood matrix and SpatialPolygonsDataFrame (region names/order)
  dat <- shape.H(shape.file)
  H <- dat$H
  map <- dat$map

  ## Attach name columns to data
  x.data <- data.frame(as.vector(x.reg.names),X)
  y.data <- data.frame(as.vector(y.reg.names),y)
  colnames(x.data)[1] <- c("NAME")
  colnames(y.data)[1] <- c("NAME")

  ## Extract the column of region names from the shapefile
  ## If no NAME column in the shapefile, the user will need to input a vector for shp.reg.names
  if(is.null(shp.reg.names)==FALSE) {
    shape.names <- shp.reg.names
  } else if ("NAME" %in% colnames(map)){
    shape.names <- map$NAME} else {
      warning("Please enter a vector of names corresponding to your shapefile")}

  ## Re-order the data to be in the same region order as the shapefile
  x.data <- x.data[order(factor(x.data$NAME,levels=shape.names)),]
  y.data <- y.data[order(factor(y.data$NAME,levels=shape.names)),]

  ## Remove region name columns for sampling
  Y <- as.vector(y.data[,-1])
  X <- as.matrix(x.data[,-1])

  ## Perform sampling using ICAR reference prior (Keefe et al., 2018)
  MCMC <- ref.MCMC(y=Y,X=X,H=H,iters,burnin,verbose,
                   tauc.start,beta.start,sigma2.start,step.tauc,
                   step.sigma2)

  ## Parameter summaries
  ref.params <- ref.summary(MCMCchain=MCMC$MCMCchain,tauc.MCMC=MCMC$tauc.MCMC,
                            sigma2.MCMC=MCMC$sigma2.MCMC,beta.MCMC=MCMC$beta.MCMC,phi.MCMC=MCMC$phi.MCMC,
                            accept.phi=MCMC$accept.phi,accept.sigma2=MCMC$accept.sigma2,
                            accept.tauc=MCMC$accept.tauc,iters,burnin)

  ## Extract parameter summaries to avoid list returns for user
  beta.median <- ref.params$beta.median
  beta.hpd <- ref.params$beta.hpd

  tauc.median <- ref.params$tauc.median
  tauc.hpd <- ref.params$tauc.hpd

  sigma2.median <- ref.params$sigma2.median
  sigma2.hpd <- ref.params$sigma2.hpd

  sigma2.accept <- ref.params$sigma2.accept
  tauc.accept <- ref.params$tauc.accept

  ## Parameter trace plots
  ref.plot(MCMC$MCMCchain,X,burnin,num.reg)

  ## Regional medians and intervals
  regions <- reg.summary(MCMC$MCMCchain,X,Y,burnin)

  fit.dist <- regions$fit.dist
  reg.medians <- regions$reg.medians
  reg.hpd <- regions$reg.hpd

  return(list(H=H,MCMC=MCMC$MCMCchain,beta.median=beta.median,
              beta.hpd=beta.hpd,tauc.median=tauc.median,tauc.hpd=tauc.hpd,
              sigma2.median=sigma2.median,sigma2.hpd=sigma2.hpd,
              tauc.accept=tauc.accept,sigma2.accept=sigma2.accept,
              fit.dist=fit.dist,reg.medians=reg.medians,reg.hpd=reg.hpd))
}

