#' Compute confidence intervals of a structural VAR(p) model using the wild bootstrap
#'
#' @param olsobj Least-squares object from running `olsvarc()`
#' @param irfobj IRF object from running `irfvar()`.
#' @param h Horizon for IRFs
#' @param nrep Positive integer for number of bootstrap replications
#' @param method String for method of computing standard errors.
#' @param bootstrap_seed Positive integer for seed
#' @param display_progress_bar Logical for whether to display the progress bar
#'
#' @return A tibble dataframe containing upper and lower ends of the confidence intervals associated with impulse response functions.
#' @import cli
#' @export
#'
#' @examples
#' y <- kilianLutkepohlCh12Figure12_5 |>
#'   dplyr::rename(oilsupply = oil_supply,
#'                 aggdemand = agg_demand) |>
#'   dplyr::select(-date)
#' var_order <- c("oilsupply", "aggdemand", "rpoil")
#' negative_shocks <- "oilsupply"
#' var_cumsum <- c("response_shock_oilsupply_oilsupply",
#'                 "response_shock_oilsupply_aggdemand",
#'                 "response_shock_oilsupply_rpoil")
#' sol <- olsvarc(y, p = 24)
#' irf <- irfvar(Ahat = sol$Ahat, B0inv = t(chol(sol$SIGMAhat)), p = sol$p, h = 15,
#'               var_order = var_order,
#'               negative_shocks = negative_shocks,
#'               var_cumsum = var_cumsum)
#' dat_ci <- bootstrap_wild(sol, irf, nrep = 500, method = "standard", standard_factor = 2.0)
#' plot(irf$irf_tidy$horizon, irf$irf_tidy$response_shock_oilsupply_oilsupply, type = "l", ylim = c(-2, 0))
#' lines(dat_ci$horizon, dat_ci$response_shock_oilsupply_oilsupply_lo, lty = "dotted", col = "blue")
#' lines(dat_ci$horizon, dat_ci$response_shock_oilsupply_oilsupply_hi, lty = "dotted", col = "blue")
bootstrap_wild <- function(olsobj, irfobj, h, nrep,
                           method = "standard", standard_factor, efron_percentile_quantiles, hall_percentile_quantiles,
                           nrep_inside_boot,
                           bootstrap_seed, display_progress_bar = TRUE) {

  # Check method is correct
  if (!(method %in% c("standard", "efron_percentile", "hall_percentile", "equal_percentile_t", "symmetric_percentile_t"))) {
    stop("Your option for constructing bootstrap confidence intervals is incorrect.")
  }

  # Deal with seed
  if (!missing(bootstrap_seed)) {
    set.seed(bootstrap_seed)
  } else {
    set.seed(676)
  }

  if (missing(h)) {
    horizon <- irfobj$h
  }

  if (missing(nrep_inside_boot)) {
    nrep_inside_boot = 500
  }

  # Rename variables for ease
  y <- olsobj$y
  tt <- nrow(y)
  pp <- olsobj$p
  KK <- olsobj$K
  TT <- tt - pp

  bootstrap_var_order <- irfobj$var_order
  bootstrap_var_cumsum <- irfobj$var_cumsum
  bootstrap_negative_shocks <-irfobj$negative_shocks

  irf_matrix <- irfobj$irfm
  virf_matrix <- matrix(irf_matrix, nrow = 1)

  IK <- diag((pp - 1) * KK)
  AP0 <- matrix(0, nrow = ((pp - 1) * KK), ncol = KK)
  A <- rbind(olsobj$Ahat, cbind(IK, AP0))

  NU <- rbind(olsobj$Vhat, matrix(0, nrow = KK * pp - KK, ncol = 1))

  U <- olsobj$Uhat # dimension K x (tt - pp), tt is the number of data and pp is the VAR length

  mIRF <- matrix(1, nrow = nrep, ncol = (KK^2) * (horizon + 1))

  # Equal-tailed percentile-t intervals
  if (method == "equal_percentile_t") {
    vIRFse <- bootstrap_wild_bootse(olsobj = olsobj, irfobj = irfobj, nrep = 500)
    # matrix of bootstrapped ts
    mbootstrapt <- matrix(NA, nrow = nrep, ncol = (KK^2) * (horizon + 1))
  }

  # Develop progress bar as needed
  if (display_progress_bar == TRUE) {
    cli::cli_progress_bar("Bootstrap reps", tot = nrep)
  }

  # Loop for boostrap replications
  for (r in 1:nrep) {

    if (display_progress_bar == TRUE) {
      cli::cli_progress_update()
    }

    rY <- matrix(NA, KK * pp, tt - pp + 1)
    rU <- matrix(0, KK * pp, tt - pp)
    rydat <- matrix(NA, tt - pp, KK)

    # Initial condition for recursively building sample
    # used to pick off vector of YY
    rpos_y0 <- floor(runif(1, min = 0, max = 1) * (tt - pp + 1)) + 1
    # Z matrix from OLS
    YY <- olsobj$Z[2:(KK * pp + 1), 1:(tt - pp), drop = FALSE]
    # y_T
    YYt <- rbind(olsobj$Y[, TT, drop = FALSE], YY[, TT, drop = FALSE])
    YY <- cbind(YY, YYt[1:(KK * pp), , drop = FALSE])

    rY0 <- YY[, rpos_y0, drop = FALSE]
    rY[1:(KK * pp), 1] <- rY0

    # No iid resampling
    # Recursive design wild bootstrap
    reta1 <- matrix(rnorm(tt - pp, mean = 0, sd = 1), nrow = 1)
    reta <- rbind(reta1, reta1, reta1)
    rU[1:KK, 1:(tt - pp)] <- U[1:KK, 1:(tt - pp)] * reta

    for (i in 2:(tt - pp + 1)) {
      rY[1:(KK * pp), i] <- NU + A %*% rY[1:(KK * pp), i - 1, drop = FALSE] + rU[1:(KK * pp), i - 1, drop = FALSE]
      rydat[i - 1, 1:KK] <- rY[1:KK, i]
    }

    # Concatenate the initial vector to rydat
    rY0_wide <- matrix(rY0, nrow = pp, ncol = KK, byrow = TRUE)
    rydat <- rbind(rY0_wide, rydat)

    rsol <- olsvarc(y = rydat, p = pp)

    if (method == "equal_percentile_t") {
      rvIRFse <- bootstrap_wild_bootse(olsobj = olsobj,
                                           irfobj = irfobj,
                                           nrep = nrep_inside_boot)
    }
    if (method == "symmetric_percentile_t") {
      rvIRFse <- bootstrap_wild_bootse(olsobj = olsobj,
                                           irfobj = irfobj,
                                           nrep = nrep_inside_boot)
    }

    rB0inv <- t(chol(rsol$SIGMAhat))
    rIRF <- irfvar(Ahat = rsol$Ahat,
                   B0inv = rB0inv,
                   p = pp,
                   h = horizon,
                   negative_shocks = bootstrap_negative_shocks,
                   var_cumsum = bootstrap_var_cumsum,
                   var_order = bootstrap_var_order)

    mIRF[r, ] <- matrix(rIRF$irfm, nrow = 1)

    if (method == "equal_percentile_t") {
      mIRF[r, ] <- ifelse(vIRFse == 0, 0, (mIRF[r, ] - virf_matrix) / vIRFse)
    } else if (method == "symmetric_percentile_t") {
      mIRF[r, ] <- ifelse(vIRFse == 0, 0, abs(mIRF[r, ] - virf_matrix) / vIRFse)
    }
  }

  if (display_progress_bar == TRUE) {
    cli::cli_progress_done()
  }

  # Compute the confidence intervals
  if (method == "standard") {
    mIRFstd <- matrix(apply(mIRF, 2, sd), nrow = KK^2, ncol = horizon + 1)
    if (missing(standard_factor)) {
      mCIlo <- irf_matrix - 2 * mIRFstd
      mCIhi <- irf_matrix + 2 * mIRFstd
    } else {
      mCIlo <- irf_matrix - standard_factor * mIRFstd
      mCIhi <- irf_matrix + standard_factor * mIRFstd
    }
  } else if (method == "efron_percentile") {
    if (missing(efron_percentile_quantiles)) {
      mIRFqtile <- apply(mIRF, 2, quantile, probs = c(0.025, 0.975))

      mCIlo <- mIRFqtile[1, ]
      mCIhi <- mIRFqtile[2, ]

      # Reshape
      mCIlo <- matrix(mCIlo, nrow = KK^2, ncol = horizon + 1)
      mCIhi <- matrix(mCIhi, nrow = KK^2, ncol = horizon + 1)
    } else {
      mIRFqtile <- apply(mIRF, 2, quantile, probs = efron_percentile_quantiles)

      mCIlo <- mIRFqtile[1, ]
      mCIhi <- mIRFqtile[2, ]

      # Reshape
      mCIlo <- matrix(mCIlo, nrow = KK^2, ncol = horizon + 1)
      mCIhi <- matrix(mCIhi, nrow = KK^2, ncol = horizon + 1)
    }
  } else if (method == "hall_percentile") {
    mIRFqtile <- apply(mIRF, 2, quantile, probs = c(0.025, 0.975))

    IRFqtile_lo <- matrix(mIRFqtile[1, ], nrow = KK^2, ncol = horizon + 1)
    IRFqtile_hi <- matrix(mIRFqtile[2, ], nrow = KK^2, ncol = horizon + 1)

    mCIlo <- 2 * irf_matrix - IRFqtile_hi
    mCIhi <- 2 * irf_matrix - IRFqtile_lo
  } else if (method == "equal_percentile_t") {
    mIRFqtile <- apply(mIRF, 2, quantile, probs = c(0.025, 0.975))

    IRFqtile_lo <- mIRFqtile[1, ]
    IRFqtile_hi <- mIRFqtile[2, ]

    mCIlo <- matrix(virf_matrix - IRFqtile_hi * vIRFse, nrow = KK^2, ncol = horizon + 1)
    mCIhi <- matrix(virf_matrix - IRFqtile_lo * vIRFse, nrow = KK^2, ncol = horizon + 1)
  } else if (method == "symmetric_percentile_t") {
    vIRFqtile <- apply(mIRF, 2, quantile, probs = c(0.95))

    mCIlo <- matrix(virf_matrix - vIRFqtile * vIRFse, nrow = KK^2, ncol = horizon + 1)
    mCIhi <- matrix(virf_matrix + vIRFqtile * vIRFse, nrow = KK^2, ncol = horizon + 1)
  }

  rownames(mCIlo) <- paste0(get_irf_names(bootstrap_var_order), "_lo")
  rownames(mCIhi) <- paste0(get_irf_names(bootstrap_var_order), "_hi")

  dat_CIlo <- as_tibble(t(mCIlo)) |>
    mutate(horizon = 0:horizon)

  dat_CIhi <- as_tibble(t(mCIhi)) |>
    mutate(horizon = 0:horizon)

  dat_CI <- dat_CIlo |>
    left_join(dat_CIhi, by = join_by(horizon))

  # Return the tibble
  return(dat_CI)
}

bootstrap_wild_bootse <- function(olsobj, irfobj, nrep) {

  if (missing(nrep)) {
    nrep <- 500
  }

  # Rename variables for ease
  y <- olsobj$y
  tt <- nrow(y)
  pp <- olsobj$p
  KK <- olsobj$K
  horizon <- irfobj$h

  bootstrap_var_order <- irfobj$var_order
  bootstrap_var_cumsum <- irfobj$var_cumsum
  bootstrap_negative_shocks <-irfobj$negative_shocks

  irf_matrix <- irfobj$irfm

  IK <- diag((pp - 1) * KK)
  AP0 <- matrix(0, nrow = ((pp - 1) * KK), ncol = KK)
  A <- rbind(olsobj$Ahat, cbind(IK, AP0))

  NU <- rbind(olsobj$Vhat, matrix(0, nrow = KK * pp - KK, ncol = 1))

  U <- olsobj$Uhat # dimension K x (tt - pp), tt is the number of data and pp is the VAR length

  mIRF <- matrix(1, nrow = nrep, ncol = (KK^2) * (horizon + 1))

  for (r in 1:nrep) {
    # print(r)

    rY <- matrix(NA, KK * pp, tt - pp + 1)
    rU <- matrix(0, KK * pp, tt - pp)
    rydat <- matrix(NA, tt - pp, KK)

    # Initial condition
    rpos_y0 <- floor(runif(1, min = 0, max = 1) * (tt - pp + 1)) + 1
    rY0 <- Y[, rpos_y0, drop = FALSE]
    rY[1:(KK * pp), 1] <- rY0

    # No iid resampling
    # Recursive design wild bootstrap
    reta1 <- matrix(rnorm(tt - pp, mean = 0, sd = 1), nrow = 1)
    reta <- rbind(reta1, reta1, reta1)
    rU[1:KK, 1:(tt - pp)] <- U[1:KK, 1:(tt - pp)] * reta

    for (i in 2:(tt - pp + 1)) {
      rY[1:(KK * pp), i] <- NU + A %*% rY[1:(KK * pp), i - 1, drop = FALSE] + rU[1:(KK * pp), i - 1, drop = FALSE]
      rydat[i - 1, 1:KK] <- rY[1:KK, i]
    }

    # Concatenate the initial vector to rydat
    rY0_wide <- matrix(rY0, nrow = pp, ncol = KK, byrow = TRUE)
    rydat <- rbind(rY0_wide, rydat)

    rsol <- olsvarc(y = rydat, p = pp)
    rB0inv <- t(chol(rsol$SIGMAhat))
    rIRF <- irfvar(Ahat = rsol$Ahat,
                   B0inv = rB0inv,
                   p = pp,
                   h = horizon,
                   negative_shocks = bootstrap_negative_shocks,
                   var_cumsum = bootstrap_var_cumsum,
                   var_order = bootstrap_var_order)

    mIRF[r, ] <- matrix(rIRF$irfm, nrow = 1)
  }

  # vIRFse <- matrix(apply(mIRF, 2, sd), nrow = KK^2, ncol = horizon + 1)
  vIRFse <- matrix(apply(mIRF, 2, sd), nrow = 1)
  return(vIRFse)
}
