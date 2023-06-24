# Fractional integrated likelihood for OLM
# @description This integrated likelihood can be found analytically so does not require
# adaptive quadrature or other integration techniques.
# 'b' is the corresponding training fraction for the fractional Bayes factor \insertCite{Porter_2023}{ref.ICAR}.
#
# @author Erica M. Porter, Christopher T. Franck, and Marco A.R. Ferreira
#
# @param X A matrix of covariates, which should include a column of 1's for models with a non-zero intercept
# @param Y A vector of responses.
# @param b Training fraction for the fractional Bayes factor (FBF) approach.
#
# @return Logarithm of the fractional integrated likelihood for the corresponding ordinary linear model (OLM).
# @export

vague.logmargin.independent.b <- function(X, Y, b=0.1) {

  p <- ncol(X)
  n <- length(Y)

  logmarginal.b <- (0.5*(n*(b-1)))*log(2*pi) + (0.5*p)*log(b) + lgamma(0.5*(n-p)) -
    lgamma(0.5*(n*b-p)) + (0.5*(p-n))*log(0.5*(t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y)) +
    (0.5*(n*b-p))*log(0.5*(b*(t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y)))

  return(logmarginal.b)
}
