#' Title
#'
#' @param olsobj
#' @param irfobj
#' @param h
#' @param nrep
#' @param bootstrap_seed
#' @param print_nrep
#'
#' @return
#' @import cli
#' @export
#'
#' @examples
bootstrap_standard <- function(olsobj, irfobj, h, nrep, bootstrap_seed = NULL, print_nrep = FALSE) {

  if (!is.null(bootstrap_seed)) {
    set.seed(bootstrap_seed)
  }

  # Renamve variables for ease
  horizon <- h
  y <- olsobj$y
  tt <- nrow(y)
  pp <- olsobj$p
  KK <- olsobj$K

  var_order <- irfobj$var_order
  var_cumsum <- irfobj$var_cumsum
  negative_shocks <-irfobj$negative_shocks
  print(negative_shocks)

  irf_matrix <- irfobj$irfm

  IK <- diag((pp - 1) * KK)
  AP0 <- matrix(0, nrow = ((pp - 1) * KK), ncol = KK)
  A <- rbind(olsobj$Ahat, cbind(IK, AP0))

  NU <- rbind(olsobj$Vhat, matrix(0, nrow = KK * pp - KK, ncol = 1))

  U <- olsobj$Uhat # dimension K x (tt - pp), tt is the number of data and pp is the VAR length

  mIRF <- matrix(1, nrow = nrep, ncol = (KK^2) * (horizon + 1))

  for (r in 1:nrep) {

    if (print_nrep == TRUE) {
      print(paste0("Bootstrap replication...", r))
    }

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
    stopifnot(nrow(rydat) == tt)

    rsol <- olsvarc(y = rydat, p = pp)
    rB0inv <- t(chol(rsol$SIGMAhat))
    rIRF <- irfvar(rsol$Ahat, rB0inv, p = pp, h = horizon, var_order = var_order, var_cumsum = var_cumsum, negative_shocks = negative_shocks)

    mIRF[r, ] <- matrix(rIRF$irfm, nrow = 1)
  }

  # Compute the confidence intervals
  mIRFstd <- matrix(apply(mIRF, 2, sd), nrow = KK^2, ncol = horizon + 1)

  mCIlo <- IRF - 2 * mIRFstd
  mCIhi <- IRF + 2 * mIRFstd

  rownames(mCIlo) <- paste0(get_irf_names(var_order_global_oil), "_lo")
  rownames(mCIhi) <- paste0(get_irf_names(var_order_global_oil), "_hi")

  dat_CIlo <- as_tibble(t(mCIlo)) |>
    mutate(horizon = 0:horizon)

  dat_CIhi <- as_tibble(t(mCIhi)) |>
    mutate(horizon = 0:horizon)

  dat_CI <- dat_CIlo |>
    left_join(dat_CIhi, by = join_by(horizon))

  # Return the tibble
  return(dat_CI)
}
