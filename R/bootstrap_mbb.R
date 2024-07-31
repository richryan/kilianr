bootstrap_mbb <- function(olsobj, irfobj, h, nrep,
                          blen,
                          method = "standard", standard_factor, efron_percentile_quantiles, hall_percentile_quantiles,
                          bootstrap_seed, display_progress_bar = TRUE) {

  # Check method is correct
  if (!(method %in% c("standard", "efron_percentile", "hall_percentile"))) {
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

  # Rename variables for ease
  y <- olsobj$y
  tt <- nrow(y)
  pp <- olsobj$p
  KK <- olsobj$K

  # Block length
  nblocks <- ceiling((tt - pp) / blen)

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
    # rU <- matrix(0, KK * pp, tt - pp)
    rU <- matrix(0, KK, nblocks * blen)
    rydat <- matrix(NA, tt - pp, KK)

    # Moving block bootstrap
    rUtilde <- matrix(0, nrow = KK, ncol = nblocks * blen)

    for (b in 0:(nblocks - 1)) {
      # draw over integers 1,...,T - blen + 1
      bupos <- ceiling(runif(1) * ((tt - pp) - blen + 1))
      print("b start")
      print(bupos)
      # indexes for allocating rUtilde
      indlo <- b * blen + 1
      indhi <- b * blen + blen
      rUtilde[, indlo:indhi] <- U[, bupos:(bupos + blen - 1)]
    }

    # Re-center
    print("start recenter")
    for (ss in 1:blen) {
      for (jj in 0:(nblocks - 1)) {
        rU[1:KK, jj * blen + ss] <- rUtilde[1:KK, jj * blen + ss] - matrix(rowMeans(U[, ss:(ss + tt - pp - blen)]), ncol = 1)
      }
    }

    print(rUtilde)
    print(rU)

    # Initial condition
    rpos_y0 <- floor(runif(1, min = 0, max = 1) * (tt - pp + 1)) + 1
    rY0 <- Y[, rpos_y0, drop = FALSE]
    rY[1:(KK * pp), 1] <- rY0

    # Moving block boostrap
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
  }

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
