# Eigenvalue Decomposition for MCMC Analysis
# @description Performs eigenvalue decomposition to return eigenvalue and eigenvector
# quantities necessary for sampling from the spatial random effects.  The quantities
# returned by this function are used as inputs in the function \code{ref.MCMC}, which implements
# the MCMC sampling algorithm proposed by Keefe et al. (2018).
#
# @importFrom Rdpack reprompt
# @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#
# @param H An adjacency matrix for the data.
# @param X A matrix of covariates.
# @param y A vector of responses.
# @param lambda A propriety parameter.
#
# @return A list of decomposition values.
#     \item{D}{Diagonal matrix of the eigenvalues of \code{H}.}
#     \item{Q}{Matrix of eigenvectors of \code{H}.}
#     \item{xi}{Vector of eigenvalues.}
#     \item{Qmat}{Matrix of eigenvectors of \code{H} with one column removed for zero eigenvalue.}
#     \item{phimat}{Covariance matrix for vector of spatial random effects.}
#     \item{Sig_phi}{Spectral decomposition of covariance depending on spatial effects.}
#     \item{X}{Matrix of covariates.}
# @export
# @examples
# ## Refer to the vignette attached to the package.
ref.input <- function(H,X,y,lambda=0) {
  if (is.matrix(X)==FALSE & is.vector(X)==TRUE) {X <- as.matrix(X, ncol=1)}

  check.mat(H)

  num.reg <- length(y)
  row.names(H) <- NULL
  Q <- eigen(H,symmetric=TRUE)$vectors
  Qmat <- Q[,1:(num.reg-1)]
  eigH <- eigen(H,symmetric=TRUE)$values
  phimat <- diag(1/sqrt(lambda + eigH[1:(num.reg-1)]))
  D <- diag(eigH)
  Sig_phi <- matrix(0,num.reg,num.reg) #initialize
  for(i in 1:(num.reg-1)){
    total <- (1/(lambda + eigH[i]))*Q[,i]%*%t(Q[,i])
    Sig_phi <- Sig_phi + total
  }
  s1 <- diag(num.reg) - X%*%solve(t(X)%*%X)%*%t(X)
  Q.star <- eigen(s1)$vectors[,1:(num.reg-ncol(X))]
  M <- t(Q.star)%*%Sig_phi%*%Q.star
  U <- eigen(M)$vectors
  L <- Q.star%*%U
  xi <- eigen(t(Q.star)%*%Sig_phi%*%Q.star)$values

  return(list(D=D,Q=Q,xi=xi,Qmat=Qmat,phimat=phimat,Sig_phi=Sig_phi,X=X))
}
