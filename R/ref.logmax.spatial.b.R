# Logarithm of the integrand over tau for the fractional integrated likelihood for ICAR model
# @description Function to help avoid numerical underflow during integration.
# @author Erica M. Porter, Christopher T. Franck, and Marco A.R. Ferreira
#
# @param X A matrix of covariates, which should include a column of 1's for models with a non-zero intercept
# @param Y A vector of responses.
# @param b Training fraction for the fractional Bayes factor (FBF) approach.
# @param H Neighborhood matrix for spatial subregions.
# @param H.spectral Spectral decomposition of neighborhood matrix, if user wants to pre-compute it to save time.
# @param Sig_phi Pseudo inverse of the neighborhood matrix, if user wants to pre-compute it to save time.
# @param logmax Value subtracted from the log fractional integrated likelihood to help avoid numerical underflow.
#
# @return Maximum log fractional integrated likelihood for the corresponding ICAR model.
# @export

ref.logmax.spatial.b <- function(tauc,b=0.1,Y,X,H,H.spectral=NULL,Sig_phi=NULL) {

  integrand <- rep(NA,length(tauc))
  n <- length(Y)
  num.reg <- length(Y)
  p <- ncol(X)
  norm.constant <- (0.5*(p-n*b))*log(2*pi) - (0.5*p)*log(b)

  # Get info from adjacency matrix H
  if(is.null(H.spectral)==TRUE){
    row.names(H) <- NULL
    H.spectral <- eigen(H,symmetric=TRUE)
  }

  Q <- H.spectral$vectors
  eigH <- H.spectral$values

  if(is.null(Sig_phi)==TRUE){
    phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
    Sig_phi <- matrix(0,num.reg,num.reg) #initialize
    for(i in 1:(num.reg-1)){
      total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
      Sig_phi <- Sig_phi + total
    }
  }

  # Projection of X and spectral decomp
  s1 <- diag(num.reg) - X%*%solve(t(X)%*%X)%*%t(X)
  s1 <- (s1+t(s1))/2
  s1.eigen <- eigen(s1, symmetric = TRUE)
  Q.star <- s1.eigen$vectors[,1:(num.reg-ncol(X))]
  M <- t(Q.star)%*%Sig_phi%*%Q.star
  M.spectral <- eigen(M, symmetric = TRUE)
  U <- M.spectral$vectors
  xi <- M.spectral$values
  L <- Q.star%*%U

  # Get chi-square quantity for approximating S2
  LtY.squared <- (t(L)%*%Y)^2

  #Y.sp <- t(Q) %*% Y
  X.sp <- t(Q) %*% X

  # Set up info for omega matrix
  for(i in 1:length(tauc)){
    aux1 <- c((1+(1/(tauc[i]*(eigH[1:(n-1)])))),1)
    omega.log.det <- sum(log(aux1))
    Ainv <- matrix(1/aux1, nrow=n, ncol=ncol(X))
    s2.tauc <- sum(LtY.squared/(1+((tauc[i]^(-1))*xi)))

    aux0 <- t(X.sp*Ainv)
    aux <- determinant(aux0 %*% X.sp,logarithm=T)$modulus[1]
    aux2 <- xi/(tauc[i]+xi)

    # Fractional integrand (with all normalizing constants)
    integrand[i] <- norm.constant - (0.5*b)*omega.log.det -
      (0.5*aux) - log(tauc[i]) +
      lgamma((0.5*(n*b-p))) + (0.5*(p-n*b))*log(0.5*(b*s2.tauc)) +
      (0.5*log((sum((aux2)^2) - (1/(n-p))*((sum((aux2)))^2))))
  }
  return(integrand)
}
