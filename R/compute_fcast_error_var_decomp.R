#' Compute how much of the prediction mean squared error is accounted for by each structural shock
#'
#' The forecast error variance decomposition of the structural model computes
#' how much of prediction mean squared error of y_{t+h} at horizon h=1,...,H is
#' accounted for by each structural shock.
#'
#' @param solvar A list that is output from running `olsvarc()`.
#' @param h A number that represents the horizon.
#' @param varpos A string that names the variable to be explained.
#' @param eps A number to check the absolute difference between percent of h-step ahead forecast error variance explained.
#'
#' @return If h is finite, a vector of numbers that represents the percent explained; otherwise, a tibble of the percent explained at different horizons.
#' @export
#'
#' @examples
#' y_stock_market <- kilianLutkepohlCh04Table4_1 |> select(all_of(var_order_stock_market))
#' sol_stock_market <- olsvarc(y_stock_market, p = 24)
#' compute_fcast_error_var_decomp(solvar = sol_stock_market, h = 3, var_name = "dd")
compute_fcast_error_var_decomp <- function(solvar, h, var_name, eps = 1e-4) {
  if (!is.infinite(h)) {
    ret <- compute_fcast_error_var_decomp_finite(solvar = solvar, h = h, var_name = var_name)
  }

  if (is.infinite(h)) {

    dat_1 <- compute_fcast_error_var_decomp_finite(solvar = solvar, h = 1, var_name = var_name)
    dat_percent_explained <- bind_cols(tibble(horizon = 1), as_tibble(dat_1))

    check <- 2 * eps + 1
    hinf <- 2
    fevd_old <- dat_1
    while (check > eps) {
      fevd_i <- compute_fcast_error_var_decomp_finite(solvar = solvar, h = hinf, var_name = var_name)
      check <- max(abs(fevd_old - fevd_i))
      print(paste0("   Starting horizon...", hinf, "...with maximum change in percentage explained = ", check))
      dat_percent_explained_i <- bind_cols(tibble(horizon = hinf), as_tibble(fevd_i))
      dat_percent_explained <- bind_rows(dat_percent_explained, dat_percent_explained_i)
      fevd_old <- fevd_i
      hinf <- hinf + 1
    }

    ret <- dat_percent_explained
  }

  return(ret)
}

compute_fcast_error_var_decomp_finite <- function(solvar, h, var_name) {
  SIGMAhat <- solvar$SIGMAhat

  varpos <- match(var_name, names(solvar$y))

  p <- solvar$p
  k <- solvar$K

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

  TH4 <- colSums(TH3)

  VC <- matrix(0, nrow = k, ncol = k)
  for (j in seq(1, k)) {
    VC[j, ] <- TH3[j, ] / TH4
  }

  ret <- t(VC[, varpos]) * 100

  if (h > 1) {
    for (i in seq(2, h, by = 1)) {
      TH <- J %*% (A %^% (i-1)) %*% t(J) %*% t(chol(SIGMAhat))
      TH <- t(TH)
      TH2 <- TH * TH
      TH3 <- TH3 + TH2
    }

    TH4 <- colSums(TH3)

    VC <- matrix(0, nrow = k, ncol = k)
    for (j in seq(1, k)) {
      VC[j, ] <- TH3[j, ] / TH4
    }

    ret <- t(VC[, varpos]) * 100
  }

  colnames(ret) <- names(solvar$y)
  return(ret)
}
