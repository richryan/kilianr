#' Compute second-stage impulse response functions
#'
#' @param y Vector or matrix of time series
#' @param x Vector or matrix of shock series
#' @param p A positive integer that represents the horizon
#' @param block_length Block length of d
#' @param standard_deviation_factor Non-negative number to scale up the standard deviation of the bootstrap draws.
#' @param nrep Number that represents the number of bootstrap draws
#' @param boot_seed Number to seed the bootstrap
#' @param cumeffect Logical to indicate whether to cummulate the impulse response
#'
#' @return A tibble with second-stage impulse response and upper- and lower-ends that define the confidence interval.
#' @import tibble
#' @import tidyr
#' @import purrr
#' @import ggplot2
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#' library(ggplot2)
#' library(lubridate)
#' y <- kilianLutkepohlCh12Figure12_5 |>
#'   dplyr::rename(oilsupply = oil_supply,
#'                 aggdemand = agg_demand) |>
#'   dplyr::select(-date)
#' var_order <- c("oilsupply", "aggdemand", "rpoil")
#' negative_shocks <- "oilsupply"
#' var_cumsum <- c("response_shock_oilsupply_oilsupply",
#'                 "response_shock_oilsupply_aggdemand",
#'                 "response_shock_oilsupply_rpoil")
#' p <- 24
#' sol <- olsvarc(y, p = p)
#' # Estimate structural shocks
#' Ehat <- solve(t(chol(sol$SIGMAhat))) %*% sol$Uhat
#' # Create tibble of structural shocks
#' dates <- kilianLutkepohlCh12Figure12_5|>
#'   slice(seq(from = p + 1, to = nrow(kilianLutkepohlCh12Figure12_5), by = 1)) |>
#'   pull(date)
#' dat_structural_shocks <- as_tibble(t(Ehat)) |>
#'   rename_with(.fn = function(x) paste0("shock_", x), .cols = everything()) |>
#'   mutate(date = dates)
#'
#' # Generate some random data to work with
#' dat_structural_shocks <- dat_structural_shocks |>
#'   mutate(somey = rnorm(n(), mean = 3, sd = 1))
#'
#' # Create 2nd-stage IRFs
#' dat_irf2 <- irf2(y = dat_structural_shocks$somey, x = dat_structural_shocks$shock_oilsupply,
#'                  p = 12, block_length = 12, standard_deviation_factor = 2, boot_seed = 1, nrep = 500)
#'
#' ggplot(data = dat_irf2) +
#'   geom_line(mapping = aes(x = horizon, y = irf2)) +
#'   geom_line(mapping = aes(x = horizon, y = irf2_lo), color = "blue", linetype = "dotted") +
#'   geom_line(mapping = aes(x = horizon, y = irf2_hi), color = "blue", linetype = "dotted")
irf2 <- function(y, x, p, block_length, standard_deviation_factor, nrep = 2000, boot_seed = 676, cumeffect = FALSE) {
  # Function for making lags
  multilag <- function(x, p) {
    lags <- 1:p
    names(lags) <- as.character(lags)
    map(lags, lag, x = x) |>
      as_tibble()
  }

  pp <- p
  datx <- as_tibble_col(x, column_name = "shock__")
  daty <- as_tibble_col(y, column_name = "y")
  dat <- bind_cols(daty, datx) |>
    mutate(across(all_of("shock__"), \(x) multilag(x, p = pp), .unpack = TRUE)) |>
    drop_na()

  y <- as.matrix(dat |> pull(y))
  X <- dat |> select(starts_with("shock__")) |> as.matrix()
  X <- cbind(matrix(1, nrow = nrow(X), ncol = 1), X)

  bhat <- solve(t(X) %*% X) %*% t(X) %*% y
  bhat_irf <- bhat[2:(pp + 1 + 1), , drop = FALSE]

  if (cumeffect == TRUE) {
    bhat_irf <- cumsum(bhat_irf)
  }

  ehat <- y - X %*% bhat

  # Block length
  ll <- block_length
  nblocks <- ceiling(length(y)/ll)

  IRF <- matrix(NA, nrow = nrep, ncol = pp + 1) # current plus lags

  set.seed(boot_seed)
  for (rr in 1:nrep) {
    ry <- matrix(NA, nrow = nblocks * ll, ncol = 1)
    rX <- matrix(NA, nrow = nblocks * ll, ncol = ncol(X))
    for (bb in 1:nblocks) {
      pos <- ceiling(runif(1) * (length(ehat) - ll))
      y_bb <- y[pos:(pos + ll - 1), , drop = FALSE]
      X_bb <- X[pos:(pos + ll - 1), , drop = FALSE]
      ry[(bb - 1) * ll + 1:ll,] <- y_bb
      rX[(bb - 1) * ll + 1:ll,] <- X_bb
    }

    rbhat <- solve(t(rX) %*% rX) %*% t(rX) %*% ry

    if (cumeffect == TRUE) {
      IRF[rr, ] <- cumsum(rbhat[2:(pp + 2)])
    } else {
      IRF[rr, ] <- rbhat[2:(pp + 2)]
    }
  }

  irfstd <- standard_deviation_factor * apply(IRF, 2, sd)

  dat_ret <- tibble(
    horizon = 0:pp,
    irf2 = as.vector(bhat_irf),
    irf2_lo = as.vector(bhat_irf - irfstd),
    irf2_hi = as.vector(bhat_irf + irfstd),
  )
  return(dat_ret)
}
