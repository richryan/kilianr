#' Compute structural impulse responses.
#'
#' @param Ahat A matrix of VAR parameters.
#' @param B0inv A matrix that corresponds to the structural impact multiplier matrix.
#' @param p A number that represents the number of lags.
#' @param h A number that represents the horizon.
#' @param var_order A string of variable names that specify the VAR ordering.
#'
#' @return A matrix of impulse responses
#' @importFrom expm %^%
#' @export
#'
#' @examples
#' IRF <- irfvar(sol$Ahat, B0inv, p, h)
irfvar <- function(Ahat, B0inv, p, h, var_order) {
  K <- nrow(B0inv)

  IK <- diag((p - 1) * K)
  AP0 <- matrix(0, nrow = ((p - 1) * K), ncol = K)
  A <- rbind(Ahat, cbind(IK, AP0))

  J <- cbind(diag(K), matrix(0, nrow = K, ncol = K * (p - 1)))

  IRF <- J %*% (A %^% 0) %*% t(J) %*% B0inv
  # To get
  # > mat_wide <- matrix(c(1,2,3,4,5,6), ncol = 3, byrow = TRUE)
  # [1, 2, 3
  #  4, 5, 6]
  #  into the matrix
  #  [1
  #   2
  #   3
  #   4
  #   5
  #   6]
  #   you need
  #   >  mat_long <- matrix(t(mat_wide), ncol = 1)
  IRF <- matrix(t(IRF), ncol = 1)

  for (i in 1:h) {
    iIRF <- J %*% (A %^% i) %*% t(J) %*% B0inv
    # See note about about using transpose in this case
    IRF <- cbind(IRF, matrix(t(iIRF), ncol = 1))
  }

  rownames(IRF) <- get_irf_names(var_order)

  dat_IRF <- as_tibble(t(IRF))
  dat_IRF <- dat_IRF |>
    mutate(horizon = 0:h) |>
    select(horizon, everything())

  return(dat_IRF)
}

get_irf_names <- function(v_names) {
  irf_names <- kronecker(v_names, v_names, FUN = paste, sep = "_")
  irf_names <- paste0("response_shock_", irf_names)
}
