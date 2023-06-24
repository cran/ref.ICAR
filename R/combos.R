
# All Combinations of a Hierarchy of Models of n Variables
# @importFrom gtools combinations
# @description Lists a matrix of combinations of 1 to n variables in ascending order.
#
# @author Chris Walsh \email{cwalsh@unimelb.edu.au}
#
# @param n an integer: the number of variables.
# @return A list containing.
#   \item{ragged}{a matrix with zeroes in empty elements and 1 in column 1, 2 in column 2 ... n in column n for full elements}
#   \item{binary}{a matrix as for ragged, but with 1 in all full elements}
# @export
combos <- function(n) {
  if (n < 2) {
    cat("\nn must be greater than 1.\n\n")
  } else if (n > 11) {
    cat("\n CAUTION! Output size increases exponentially with n. \n")
    cat(" Result for n = ", n, "will be a ", 2^n - 1, "*", n, " matrix.\n")
    cat(" Do you really want to proceed?\n")
    choice <- utils::menu(c("Proceed", "Quit"))
    if (choice == 1)
      combos1(n) else {
      }
  } else combos1(n)
}

combos1 <- function(n) {
  x <- cbind(gtools::combinations(n, 1, 1:n), array(0, dim = c(n, n - 1)))
  for (i in 2:n) {
    nc <- dim(gtools::combinations(n, i, 1:n))[1]
    x <- rbind(x, cbind(gtools::combinations(n, i, 1:n),
                        array(0, dim = c(nc, n - i))))
  }
  len <- dim(x)[1]
  x.index <- cbind(as.vector(1:len), as.vector(x))
  x.index <- cbind(x.index[, 1][x.index[, 2] > 0],
                   x.index[, 2][x.index[, 2] > 0])
  x.bin <- array(0, dim = c(len, n))
  x.bin[x.index] <- 1
  list(ragged = x, binary = x.bin)
}
