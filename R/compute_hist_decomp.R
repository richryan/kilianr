#' Compute historical decomposition
#'
#' @param olsvarc_obj Object created by `olsvarc()`
#' @param rownames Names for returned T - p tibble (most likely `date` from the original data)
#' @param var_order Order of variables for structural VAR model
#'
#' @return
#' @export
#'
#' @examples
compute_hist_decomp <- function(olsvarc_obj, rownames, var_order) {

  # Rename for ease
  sol <- olsvarc_obj

  B0inv <- t(chol(ols$SIGMAhat))
  B0 <- solve(B0inv)

  tt <- nrow(ols$y)
  pp <- sol$p

  irf <- irfvar(
    Ahat = sol$Ahat,
    B0inv = B0inv,
    var_order = var_order,
    p = sol$p,
    h = tt - pp - 1
  )

  Ehat <- B0 %*% sol$Uhat
  IRF <- irf$irfm

  hdecomp_names <-
    str_replace(hdecomp_names_, "response_shock", "hd_series_shock")

  hdecomp <- matrix(0, nrow = nrow(IRF), ncol = ncol(IRF))
  rownames(hdecomp) <- hdecomp_names

  for (series in var_order) {
    for (shock in var_order) {
      series_shock <- paste(series, shock, sep = "_")
      hd_name <- paste0("hd_series_shock_", series_shock)
      irf_name <- paste0("response_shock_", series_shock)

      series_shock_yhat <- vector(mode = "double", length = tt - pp)

      for (i in 1:(tt - pp)) {
        series_shock_yhat[i] <-
          IRF[irf_name, 1:i, drop = FALSE] %*% t(Ehat[shock, seq(i, 1, by = -1), drop = FALSE])
      }

      # print(series_shock_yhat)
      hdecomp[hd_name, ] <- series_shock_yhat
    }
  }

  return(hdecomp)
}

