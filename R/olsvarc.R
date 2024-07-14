#' Estimate a VAR model with an intercept using OLS.
#'
#' Details are provided on page 31 of Kilian, Lutz, and Helmut Lütkepohl. *Structural Vector Autoregressive Analysis.*
#'
#' @param y A matrix of data.
#' @param p A positive integer, indicating the number of lags.
#'
#' @returns A list with components\tabular{ll}{
#'    \code{Vhat} \tab A numeric vector of parameter estimates \cr
#'    \tab \cr
#'    \code{Ahat} \tab OLS estimates of A_i, i = 1,..., p, where each A_i is K by K \cr
#'    \tab \cr
#'    \code{SIGMAhat} \tab Consistent estimator of the innovation covariance matrix \cr
#'    \tab \cr
#'    \code{Uhat} \tab Estimated least-squares residuals \cr
#'    \tab \cr
#'    \code{Z} \tab Constructed data matrix \cr
#' }
#' @export
#'
#' @examples
#' sol <- olsvarc(y, p)
olsvarc <- function(y, p) {
  pp <- p
  tt <- nrow(y)
  TT <- tt - pp
  KK <- ncol(y)

  X <- t(as.matrix(y))
  Y <- X[, (pp + 1):tt]
  # Create Z matrix by row-binding a matrix below,
  # starting with vector of ones
  Z <- matrix(1, nrow = 1, ncol = TT)
  for (i in 1:pp) {
    Z <- rbind(Z, X[, (pp + 1 - i):(tt - i)])
  }

  Ahat <- Y %*% t(Z) %*% solve(Z %*% t(Z))

  Uhat <- Y - Ahat %*% Z
  SIGMAhat <- Uhat %*% t(Uhat) / (TT - KK * pp - 1)

  Vhat <- Ahat[, 1, drop = FALSE]
  Ahat <- Ahat[, 2:(KK * pp + 1)]

  return(list(Vhat=Vhat, Ahat=Ahat, SIGMAhat = SIGMAhat, Uhat=Uhat, Z=Z, Y=Y, p=pp))
}
