# Calculating the Moore-Penrose Inverse of a Matrix
# @description Calculates the Moore-Penrose inverse.
# @author Marco Ferreira
# @param Sigma The matrix to invert.
#
# @return The Moore-Penrose inverse
#     \item{Sigma.M.P.Inverse}{Moore-Penrose inverse of Sigma}
# @export
M.P.Inverse <- function(Sigma)
{
  # M.A.R. Ferreira (2009)
  # The spectral decomposition of Sigma = P D P', eliminating the columns corresponding to null eigenvalues is:
  Sigma.D <- eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values
  number.redundant <- sum((Sigma.D)^2 < 0.00000000001)
  dim.sig <- length(Sigma.D) - number.redundant
  Sigma.D <- Sigma.D[1:dim.sig]
  Sigma.P <- svd(Sigma,nu=dim.sig,nv=dim.sig)$u
  # The Moore-Penrose inverse of Sigma is:
  Sigma.M.P.Inverse <- Sigma.P %*% diag(1.0/Sigma.D) %*% t(Sigma.P)
  Sigma.M.P.Inverse
}
