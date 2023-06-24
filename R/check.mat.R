# Checking a Neighborhood Matrix for Areal Data
# @description Checks that the Neighborhood matrix, \code{H}, is correct for analysis using the
# objective ICAR model \insertCite{Keefe_2018}{ref.ICAR}.  Specifically, \code{H} must be symmetric and contain
# only integer-valued elements.  The diagonals correspond to the number of neighbors for each region,
# and each off-diagonal element equals -1 if the corresponding regions are neighbors, and 0 otherwise.
#
# @author Erica M. Porter, Matthew J. Keefe, Christopher T. Franck, and Marco A.R. Ferreira
#
# @param H neighborhood matrix for the data regions.
#
# @return Warnings on the neighborhood matrix, if any.
# @export
# @examples
# #### The following should give warnings on H.
#
# H <- matrix(c(1,-1,0,0,-1,3,-1,0,0,-1,1,0,0,-1,0,1),nrow=4,ncol=4)
#
# \donttest{check.mat(H)}
#
# H <- matrix(c(1,-1,0,0,-1,2.4,-1,0,0,-1,1.2,0,0,-1,0,1),nrow=4,ncol=4)
#
# \donttest{check.mat(H)}

check.mat <- function(H) {
  if (isSymmetric(H) == FALSE) {stop("H must be a symmetric matrix")}
  if (all(diag(H)>0) == FALSE) {stop("The diagonal elements of H must be positive")}
  if (all(diag(H) == floor(diag(H))) == FALSE) {stop("The diagonal elements of H must be integers")}
  if (mean(diag(H) == -apply(H-diag(diag(H)),1,sum)) != 1 | mean(diag(H) == -apply(H-diag(diag(H)),2,sum)) != 1) {stop("The diagonal elements of H should equal the number of neighbors for each region")}
  if (sum(eigen(H)$values < 1e-5) > 1) {stop("The specified region must be contiguous for this analysis")}
}
