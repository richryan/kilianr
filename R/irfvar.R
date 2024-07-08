#' Compute structural impulse responses.
#'
#' @param Ahat A matrix of VAR parameters.
#' @param B0inv A matrix that corresponds to the structural impact multiplier matrix.
#' @param p A number that represents the number of lags.
#' @param h A number that represents the horizon.
#'
#' @return A matrix of impulse responses
#' @importFrom expm %^%
#' @export
#'
#' @examples
#' IRF <- irfvar(sol$Ahat, B0inv, p, h)
irfvar <- function(Ahat, B0inv, p, h) {
  K <- nrow(B0inv)

  IK <- diag((p - 1) * K)
  AP0 <- matrix(0, nrow = ((p - 1) * K), ncol = K)
  A <- rbind(Ahat, cbind(IK, AP0))

  J <- cbind(diag(K), matrix(0, nrow = K, ncol = K * (p - 1)))

  IRF <- J %*% (A %^% 0) %*% t(J) %*% B0inv
  # To get
  # > mat_wide <- matrix(c(1,2,3,4,5,6), ncol = 3, byrow = TRUE)
  # [1, 2, 3
  #  4, 5, 6]
  #  into the matrix
  #  [1
  #   2
  #   3
  #   4
  #   5
  #   6]
  #   you need
  #   >  mat_long <- matrix(t(mat_wide), ncol = 1)
  IRF <- matrix(t(IRF), ncol = 1)

  for (i in 1:h) {
    iIRF <- J %*% (A %^% i) %*% t(J) %*% B0inv
    # See note about about using transpose in this case
    IRF <- cbind(IRF, matrix(t(iIRF), ncol = 1))
  }

  return(IRF)
}
