#' MCMC for Reference Prior on an Intrinsic Conditional Autoregressive Random Effects Model for Areal Data
#'
#' @description Implements the Metropolis-within-Gibbs sampling algorithm proposed by Keefe et al. (2018),
#' to perform posterior inference for the intrinsic conditional autoregressive model with spatial
#' random effects.
#'
#' @references
#' \insertRef{keefe2018}{ref.ICAR}
#'
#' \insertRef{Keefe_2018}{ref.ICAR}
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

    num.reg<-length(y)

    #Eigenvalue decomposition for sampling
    if (is.matrix(X)==FALSE & is.vector(X)==TRUE) {X <- as.matrix(X, ncol=1)}

    check.mat(H)
    row.names(H) <- NULL
    Q <- eigen(H,symmetric=TRUE)$vectors
    Qmat <- Q[,1:(num.reg-1)]
    eigH <- eigen(H,symmetric=TRUE)$values
    phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
    D <- diag(eigH)
    Sig_phi <- matrix(0,num.reg,num.reg) #initialize
    for(i in 1:(num.reg-1)){
        total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
        Sig_phi <- Sig_phi + total
    }
    s1 <- diag(num.reg) - X%*%solve(t(X)%*%X)%*%t(X)
    Q.star <- eigen(s1)$vectors[,1:(num.reg-ncol(X))]
    M <- t(Q.star)%*%Sig_phi%*%Q.star
    U <- eigen(M)$vectors
    L <- Q.star%*%U
    xi <- eigen(t(Q.star)%*%Sig_phi%*%Q.star)$values

    #Initialize Metropolis Algorithm
    tauc.MCMC<-matrix(0,nrow=iters,ncol=1)
    sigma2.MCMC<-matrix(0,nrow=iters,ncol=1)
    beta.MCMC<-matrix(0,nrow=iters,ncol=ncol(X))
    phi.MCMC<-matrix(0,nrow=iters,ncol=num.reg)
    tauc.MCMC[1,]<-tauc.start
    sigma2.MCMC[1,]<-sigma2.start
    phi.MCMC[1,]<-scale((1:num.reg)/100,center=TRUE,scale=FALSE)
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

        N0<-(-1/2)*(num.reg*log(prop.sigma2) + sum(log((1+(1/(prop.tauc*(diag(D)[1:(num.reg-1)]))))))) -
            ((1/(2*prop.sigma2))*t(y-X%*%curr.beta)%*%(Q%*%diag(c((1+(1/(prop.tauc*(diag(D)[1:(num.reg-1)]))))^(-1),1))%*%t(Q))%*%(y-X%*%curr.beta)) -
            log(prop.sigma2) + ((0.5)*log(((sum((xi/(prop.tauc+xi))^2)) - ((1/(num.reg - ncol(X)))*(sum((xi/(prop.tauc+xi))))^2)))) + log(prop.tauc) +log(prop.sigma2) - log(prop.tauc)
        D0<-(-1/2)*(num.reg*log(curr.sigma2) + sum(log((1+(1/(curr.tauc*(diag(D)[1:(num.reg-1)]))))))) -
            ((1/(2*curr.sigma2))*t(y-X%*%curr.beta)%*%(Q%*%diag(c((1+(1/(curr.tauc*(diag(D)[1:(num.reg-1)]))))^(-1),1))%*%t(Q))%*%(y-X%*%curr.beta)) -
            log(curr.sigma2) + ((0.5)*log(((sum((xi/(curr.tauc+xi))^2)) - ((1/(num.reg - ncol(X)))*(sum((xi/(curr.tauc+xi))))^2)))) + log(curr.tauc) + log(curr.sigma2) - log(curr.tauc)

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
        Ainv<-(Q%*%diag(c((1+(1/(curr.tauc*(diag(D)[1:(num.reg-1)]))))^(-1),1))%*%t(Q))
        mu<-solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%y
        V<-curr.sigma2*solve(t(X)%*%Ainv%*%X)
        beta.MCMC[i,]<-t(rmvnorm(1,mean=mu,sigma=V))
        curr.beta<-beta.MCMC[i,]
        accept.beta<-accept.beta + 1

        #Sample Phis from full conditional
        mu1<-X%*%curr.beta
        mu2<-rep(0,num.reg)
        Sig_11<-curr.sigma2*(diag(num.reg) + ((1/curr.tauc)*Sig_phi))
        Sig_12<-(curr.sigma2/curr.tauc)*Sig_phi
        Sig_21<-Sig_12
        Sig_22<-Sig_12
        phi.fc.mean<-mu2 + Sig_21%*%M.P.Inverse(Sig_11)%*%(y - mu1)
        phi.fc.var <- Sig_22 - Sig_21%*%M.P.Inverse(Sig_11)%*%Sig_12
        #make variance matrix symmetric to account for machine precision
        phi.fc.var.sym <- (phi.fc.var + t(phi.fc.var))/2

        prop.phi<-t(rmvnorm(1,mean=phi.fc.mean,sigma=phi.fc.var.sym,method="svd"))

        phi.MCMC[i,]<-prop.phi
        accept.phi<-accept.phi+1

        if(verbose==TRUE) {
        if((100*(i/iters))%%5==0){print(paste(100*(i/iters),'% complete',' at ',date(),sep=''))}
        }
    }

    MCMCchain <- cbind(tauc.MCMC,sigma2.MCMC,beta.MCMC,phi.MCMC)

    return(list(MCMCchain=MCMCchain,tauc.MCMC=tauc.MCMC,sigma2.MCMC=sigma2.MCMC,beta.MCMC=beta.MCMC,
                phi.MCMC=phi.MCMC,accept.phi=accept.phi,accept.sigma2=accept.sigma2,accept.tauc=accept.tauc))
}
