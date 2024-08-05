#' Compute confidence intervals using the standard bootstrap.
#'
#' @param olsobj An object created by running `olsvarc()`
#' @param irfobj An object created by running `irfvar()`
#' @param h A positive integer indicating the horizon of the impulse response
#' @param nrep Number of boostrap replications
#' @param method String for method of computing standard errors
#' @param standard_factor A number for scaling the estimated standard deviations
#' @param efront_percentile_quantiles A number
#' @param hall_percentile_quantiles A list of quantiles
#' @param nrep_inside_boot Number of boostrap replications for the inside
#' @param bootstrap_seed A number for setting the seed for the inside boostrap
#' @param display_progress_bar Logical for indicating whether the progress bar should be displayed
#'
#' @return A tibble of confidence intervals
#' @import cli
#' @export
#'
#' @examples
#' y_global_oil <- kilianLutkepohlCh12Figure12_5 |>
#'   mutate(across(c(oil_supply, agg_demand, rpoil), \(x) detrendcl(x, tt = "constant"))) |>
#'   rename(oilsupply = oil_supply,
#'          aggdemand = agg_demand) |>
#'   select(-date) |>
#'   as.matrix()
#'
#' sol_global_oil <- olsvarc(y_global_oil, p = 24)
#' var_order_global_oil <- c("oilsupply", "aggdemand", "rpoil")
#' var_cumsum_global_oil <- c("response_shock_oilsupply_oilsupply",
#'                            "response_shock_oilsupply_aggdemand",
#'                            "response_shock_oilsupply_rpoil")
#' negative_shocks_global_oil <- c("oilsupply")
#'
#' irf_global_oil <- irfvar(Ahat = sol_global_oil$Ahat,
#'                          B0inv = t(chol(sol_global_oil$SIGMAhat)),
#'                          p = sol_global_oil$p, h = 15,
#'                          negative_shocks = negative_shocks_global_oil,
#'                          var_cumsum = var_cumsum_global_oil,
#'                          var_order = var_order_global_oil)
#' dat_irf_global_oil <- irf_global_oil$irf_tidy
bootstrap_standard <- function(olsobj, irfobj, h, nrep,
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
    vIRFse <- bootstrap_standard_bootse(olsobj = olsobj, irfobj = irfobj, nrep = 500)
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

    # iid resampling
    riid_index <- floor(runif(tt - pp, min = 0, max = 1) * (tt - pp)) + 1
    # print(riid_index)
    rU[1:KK, 1:(tt - pp)] <- U[, riid_index]

    # Equation (2.2.6) in Kilian and Lutkepohl (2017)
    for (i in 2:(tt - pp + 1)) {
      rY[1:(KK * pp), i] <- NU + A %*% rY[1:(KK * pp), i - 1, drop = FALSE] + rU[1:(KK * pp), i - 1, drop = FALSE]
      rydat[i - 1, 1:KK] <- rY[1:KK, i]
    }

    # Concatenate the initial vector to rydat
    rY0_wide <- matrix(rY0, nrow = pp, ncol = KK, byrow = TRUE)
    rydat <- rbind(rY0_wide, rydat)

    rsol <- olsvarc(y = rydat, p = pp)

    if (method == "equal_percentile_t") {
      rvIRFse <- bootstrap_standard_bootse(olsobj = olsobj,
                                           irfobj = irfobj,
                                           nrep = nrep_inside_boot)
    }
    if (method == "symmetric_percentile_t") {
      rvIRFse <- bootstrap_standard_bootse(olsobj = olsobj,
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

bootstrap_standard_bootse <- function(olsobj, irfobj, nrep) {

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

    # iid resampling
    riid_index <- floor(runif(tt - pp, min = 0, max = 1) * (tt - pp)) + 1
    # print(riid_index)
    rU[1:KK, 1:(tt - pp)] <- U[, riid_index]

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
