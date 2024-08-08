#' Compute historical decomposition
#'
#' @param olsvarc_obj Object created by `olsvarc()`
#' @param date The `date` from the original data (trimmed to T - p)
#' @param var_order Order of variables for structural VAR model
#'
#' @return A tibble dataframe that shows how
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom stringr str_replace
#' @export
#'
#' @examples
#' y <- kilianLutkepohlCh12Figure12_5 |>
#'   dplyr::rename(oilsupply = oil_supply,
#'                 aggdemand = agg_demand) |>
#'   dplyr::select(-date)
#' sol <- olsvarc(y, p = 24)
#' hd <- compute_hist_decomp(sol, date = kilianLutkepohlCh12Figure12_5$date, var_order = c("oilsupply", "aggdemand", "rpoil"))
#' plot(hd$date, hd$hd_series_shock_oilsupply_oilsupply, type = "l")
#' plot(hd$date, hd$hd_series_shock_aggdemand_oilsupply + hd$hd_series_shock_aggdemand_aggdemand + hd$hd_series_shock_aggdemand_rpoil, type = "l", col = "red")
#' points(kilianLutkepohlCh12Figure12_5$date, kilianLutkepohlCh12Figure12_5$agg_demand, type = "l")
#'
#' # Another example ---------------------------------------------------------
#' y <- kilianLutkepohlCh02Macro |> select(drgnp, irate, infl)
#' m <- olsvarc(y, p = 4)
#' hd <- compute_hist_decomp(m, date = kilianLutkepohlCh02Macro$date, var_order = c("drgnp", "irate", "infl"))
compute_hist_decomp <- function(olsvarc_obj, date, var_order) {

  # Rename for ease
  sol <- olsvarc_obj

  B0inv <- t(chol(sol$SIGMAhat))
  # Inverse of B0inv = B0
  B0 <- solve(B0inv)

  tt <- nrow(sol$y)
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


  hdecomp_names_ <- get_irf_names(var_order)

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

  # Return tibble
  hdecomp_t <- t(hdecomp)
  dat_hdecomp <- as_tibble(hdecomp_t) |>
    mutate(date = date[(pp+1):tt])

  return(dat_hdecomp)
}

