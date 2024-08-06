#' Compute how much of the prediction mean squared error is accounted for by each structural shock
#'
#' The forecast error variance decomposition of the structural model computes
#' how much of prediction mean squared error of y_{t+h} at horizon h=1,...,H is
#' accounted for by each structural shock.
#'
#' @param solvar A list that is output from running `olsvarc()`.
#' @param k A number that represents the number of series: y_t is k by 1.
#' @param p A number that represents the order of the VAR process; that is, the number of lags.
#' @param h A number that represents the horizon.
#' @param varpos A number that represents position of the series.
#'
#' @return A number that represents the percent
#' @export
#'
#' @examples
#' myf1 <- fun_compute_fcast_error_decomp(solvar = solvar, k = k, p = p, h = 1, varpos = 4)
compute_fcast_error_decomp <- function(solvar, k, p, h, varpos) {
  SIGMAhat <- solvar$SIGMAhat

  Ahat <- solvar$Ahat

  Ik <- diag((p - 1) * k)
  AP0 <- matrix(0, nrow = ((p - 1) * k), ncol = k)
  A <- rbind(Ahat, cbind(Ik, AP0))

  J1 <- diag(1, nrow = k, ncol = k)
  J2 <- matrix(0, nrow = k, ncol = k * (p - 1))
  J <- cbind(J1, J2)

  TH1 <- J %*% (A %^% 0) %*% t(J)
  TH <- TH1 %*% t(chol(SIGMAhat))
  TH <- t(TH)
  TH2 <- TH * TH
  TH3 <- TH2

  if (h > 1) {
    for (i in seq(2, h, by = 1)) {
      TH <- J %*% (A %^% (i-1)) %*% t(J) %*% t(chol(SIGMAhat))
      TH <- t(TH)
      TH2 <- TH * TH
      TH3 <- TH3 + TH2
    }
  }

  TH4 <- colSums(TH3)

  VC <- matrix(0, nrow = k, ncol = k)
  for (j in seq(1, k)) {
    VC[j, ] <- TH3[j, ] / TH4
  }

  return(t(VC[, k]) * 100)
}
