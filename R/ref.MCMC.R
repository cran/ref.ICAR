#' MCMC for Reference Prior on an Intrinsic Conditional Autoregressive Random Effects Model for Areal Data
#'
#' @description Implements the Metropolis-within-Gibbs sampling algorithm proposed by Ferreira et al. (2021),
#' to perform posterior inference for the intrinsic conditional autoregressive model with spatial
#' random effects.  This algorithm uses the spectral domain for the hierarchical model to create the
#' Spectral Gibbs Sampler (SGS), which provides notable speedups to the MCMC algorithm proposed by Keefe et al (2019).
#'
#' @references
#' \insertRef{keefe2018}{ref.ICAR}
#'
#' \insertRef{Keefe_2018}{ref.ICAR}
#'
#' \insertRef{Ferreira_2021}{ref.ICAR}
#'
#' @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @importFrom stats rnorm runif
#' @importFrom mvtnorm rmvnorm
#' @param y Vector of responses.
#' @param X Matrix of covariates.  This should include a column of 1's for models with a non-zero intercept.
#' @param H The neighborhood matrix.
#' @param iters Number of MCMC iterations to perform.  Defaults to 10,000.
#' @param burnin Number of MCMC iterations to discard as burn-in.  Defaults to 5,000.
#' @param verbose If FALSE, MCMC progress is not printed.
#' @param tauc.start Starting value for the spatial dependence parameter.
#' @param beta.start Starting value for the vector of fixed effect regression coefficients.
#' @param sigma2.start Starting value for the variance of the unstructured random effects.
#' @param step.tauc Step size for the spatial dependence parameter
#' @param step.sigma2 Step size for the variance of the unstructured random effects.
#'
#' @return A list containing MCMC chains and parameter summaries.
#'     \item{MCMCchain}{Matrix of MCMC chains.}
#'     \item{tauc.MCMC}{MCMC chains for the spatial dependence parameter.}
#'     \item{sigma2.MCMC}{MCMC chains for the variance of the unstructured random effects.}
#'     \item{phi.MCMC}{MCMC chains for the spatial random effects.}
#'     \item{beta.MCMC}{MCMC chains for the fixed effect regression coefficients.}
#'     \item{accept.sigma2}{Final acceptance number for variance of the unstructured random effects.}
#'     \item{accept.tauc}{Final acceptance number for spatial dependence parameter.}
#'     \item{accept.phi}{Final acceptance number for spatial random effects.}
#' @export
#' @examples
#' #### Fit the model for simulated areal data on a grid ####
#'
#' ### Load extra libraries
#' library(sp)
#' library(methods)
#' library(spdep)
#' library(mvtnorm)
#'
#' ### Generate areal data on a grid
#' rows=5; cols=5
#' tauc=1
#' sigma2=2; beta=c(1,5)
#'
#' ### Create grid
#' grid <- GridTopology(c(1,1), c(1,1), c(cols,rows))
#' polys <- as(grid, "SpatialPolygons")
#' spgrid <- SpatialPolygonsDataFrame(polys,data=data.frame(row.names=row.names(polys)))
#'
#' ### Create neighborhood matrix
#' grid.nb <- poly2nb(spgrid,queen=FALSE)
#' W <- nb2mat(grid.nb, style="B")
#'
#' ### Put spatially correlated data in grid
#' p <- length(beta)
#' num.reg <- (rows*cols)
#' if(p>1){x1<-rmvnorm(n=num.reg,mean=rep(0,p-1),sigma=diag(p-1))} else{x1<-NULL}
#' X <- cbind(rep(1,num.reg),x1)
#' Dmat <- diag(apply(W,1,sum))
#' H <- Dmat - W
#' row.names(H) <- NULL
#'
#' ### Obtain true response vector
#' theta_true <- rnorm(num.reg,mean=0,sd=sqrt(sigma2))
#' Q <- eigen(H,symmetric=TRUE)$vectors
#' eigH <- eigen(H,symmetric=TRUE)$values
#' D <- diag(eigH)
#' Qmat <- Q[,1:(num.reg-1)]
#' phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
#' z <- t(rmvnorm(1,mean=rep(0,num.reg-1),sigma=diag(num.reg-1)))
#' phi_true <- sqrt((1/tauc)*sigma2)*(Qmat%*%phimat%*%z)
#' Y <- X%*%beta + theta_true + phi_true
#'
#' ### Fit the model
#' set.seed(5432)
#' \donttest{model <- ref.MCMC(y=Y,X=X,H=H,iters=15000,burnin=5000,verbose=TRUE,tauc.start=.1,beta.start=-1,
#' sigma2.start=.1,step.tauc=0.5,
#'       step.sigma2=0.5)
#'       }
#'
#' #### Small example for checking
#' model <- ref.MCMC(y=Y,X=X,H=H,iters=1000,burnin=50,verbose=TRUE,tauc.start=.1,beta.start=-1,
#' sigma2.start=.1,step.tauc=0.5,
#'       step.sigma2=0.5)
ref.MCMC <- function(y,X,H,iters=10000,burnin=5000,verbose=TRUE,
                     tauc.start=1,beta.start=1,sigma2.start=1,step.tauc=0.5,
                     step.sigma2=0.5){

  if (is.matrix(X)==FALSE & is.vector(X)==TRUE) {X <- as.matrix(X, ncol=1)}

  sample.size<-length(y)

  #Eigenvalue decomposition for sampling
  if (is.matrix(X)==FALSE & is.vector(X)==TRUE) {X <- as.matrix(X, ncol=1)}

  #    check.mat(H)
  row.names(H) <- NULL
  H.spectral <- eigen(H,symmetric=TRUE)
  Q <- H.spectral$vectors
  Qmat <- Q[,1:(sample.size-1)]
  eigH <- H.spectral$values
  Sig_phi <- matrix(0,sample.size,sample.size) #initialize
  for(i in 1:(sample.size-1)){
    total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
    Sig_phi <- Sig_phi + total
  }
  s1 <- diag(sample.size) - X%*%solve(t(X)%*%X)%*%t(X)
  Q.star <- eigen(s1)$vectors[,1:(sample.size-ncol(X))]
  M <- t(Q.star)%*%Sig_phi%*%Q.star
  U <- eigen(M)$vectors
  L <- Q.star%*%U
  xi <- eigen(t(Q.star)%*%Sig_phi%*%Q.star)$values


  # Spectral transformation
  y.sp <- t(Q) %*% y
  X.sp <- t(Q) %*% X

  #Initialize Metropolis Algorithm
  tauc.MCMC<-matrix(0,nrow=iters,ncol=1)
  sigma2.MCMC<-matrix(0,nrow=iters,ncol=1)
  beta.MCMC<-matrix(0,nrow=iters,ncol=ncol(X))
  phi.MCMC<-matrix(0,nrow=iters,ncol=sample.size)
  tauc.MCMC[1,]<-tauc.start
  sigma2.MCMC[1,]<-sigma2.start
  phi.MCMC[1,]<-scale((1:sample.size)/100,center=TRUE,scale=FALSE)
  beta.MCMC[1,]<-beta.start

  #initialize acceptance rates
  accept.phi<-1
  accept.tauc<-1
  accept.sigma2<-1
  accept.beta<-1

  for(i in 2:iters){

    curr.sigma2<-sigma2.MCMC[i-1,]
    curr.tauc<-tauc.MCMC[i-1,]
    curr.phi<-phi.MCMC[i-1,]
    curr.beta<-beta.MCMC[i-1,]

    #Metropolis step for sigma2 and tauc jointly
    logprop.sigma2<-rnorm(1,mean=log(curr.sigma2),sd=step.sigma2) #Propose log(sigma2) using normal distribution
    prop.sigma2<-exp(logprop.sigma2)
    logprop.tauc<-rnorm(1,mean=log(curr.tauc),sd=step.tauc) #Propose log(tauc) using normal distribution
    prop.tauc<-exp(logprop.tauc)

    aux <- y.sp - X.sp%*%curr.beta

    N0<-(-1/2)*(sample.size*log(prop.sigma2) + sum(log((1+(1/(prop.tauc*(eigH[1:(sample.size-1)]))))))) -
      (1/(2*prop.sigma2))*(t(aux * c((1+(1/(prop.tauc*eigH[1:(sample.size-1)])))^(-1),1))%*% aux) -
      log(prop.sigma2) + ((0.5)*log(((sum((xi/(prop.tauc+xi))^2)) - ((1/(sample.size - ncol(X)))*(sum((xi/(prop.tauc+xi))))^2)))) + log(prop.tauc) +log(prop.sigma2) - log(prop.tauc)

    D0<-(-1/2)*(sample.size*log(curr.sigma2) + sum(log((1+(1/(curr.tauc*(eigH[1:(sample.size-1)]))))))) -
      (1/(2*curr.sigma2))*(t(aux *c((1+(1/(curr.tauc*eigH[1:(sample.size-1)])))^(-1),1))%*% aux) -
      log(curr.sigma2) + ((0.5)*log(((sum((xi/(curr.tauc+xi))^2)) - ((1/(sample.size - ncol(X)))*(sum((xi/(curr.tauc+xi))))^2)))) + log(curr.tauc) + log(curr.sigma2) - log(curr.tauc)

    if(log(runif(1,0,1))<(N0-D0)){
      sigma2.MCMC[i,]<-prop.sigma2
      accept.sigma2<-accept.sigma2 + 1
      tauc.MCMC[i,]<-prop.tauc
      accept.tauc<-accept.tauc + 1
    }
    else{
      sigma2.MCMC[i,]<-curr.sigma2
      tauc.MCMC[i,]<-curr.tauc
    }
    curr.sigma2<-sigma2.MCMC[i,]
    curr.tauc<-tauc.MCMC[i,]

    #Gibbs step for betas
    Ainv <- matrix(c((1+(1/(curr.tauc*(eigH[1:(sample.size-1)]))))^(-1),1), nrow=sample.size, ncol=ncol(X))
    aux0 <- t(X.sp*Ainv)
    aux <- solve(aux0 %*% X.sp)
    mu <- aux %*% aux0 %*% y.sp
    V<-curr.sigma2 * aux
    beta.MCMC[i,]<-t(rmvnorm(1,mean=mu,sigma=V))
    curr.beta<-beta.MCMC[i,]
    accept.beta<-accept.beta + 1

    # Simulate spectral random effects
    #        curr.xi <- c(rnorm(n = (sample.size-1),
    #                           mean = (y.sp[1:(sample.size-1)] - X.sp[1:(sample.size-1),] %*%
    #                                       curr.beta)/(1+curr.tauc*eigH[1:(sample.size-1)]),
    #                            sd = sqrt(1/(1+curr.tauc*eigH[1:(sample.size-1)]))
    #                           )
    #                      , 0)
    eigH[sample.size] = 0.0
    curr.xi <- rnorm(n = sample.size,
                     mean = (y.sp - X.sp %*% curr.beta)/(1+curr.tauc*eigH),
                     sd = sqrt(1/(1+curr.tauc*eigH)))

    # Simulate spatial random effects
    phi.MCMC[i,] <- Q %*% curr.xi

    if(verbose==TRUE) {
      if((100*(i/iters))%%5==0){print(paste(100*(i/iters),'% complete',' at ',date(),sep=''))}
    }
  }

  MCMCchain <- cbind(tauc.MCMC,sigma2.MCMC,beta.MCMC,phi.MCMC)

  return(list(MCMCchain=MCMCchain,tauc.MCMC=tauc.MCMC,sigma2.MCMC=sigma2.MCMC,beta.MCMC=beta.MCMC,
              phi.MCMC=phi.MCMC,accept.phi=accept.phi,accept.sigma2=accept.sigma2,accept.tauc=accept.tauc))
}
