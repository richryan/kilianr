#' Remove the mean (constant) or linear trend (linear) from a vector.
#'
#' See the pracma package for more robust commands.
#'
#' @param x A vector considered a time series.
#' @param tt A string that indicates trend type, "constant" or "linear".
#'
#' @return A vector that has removed the mean or linear trend from `x`. That is,
#' for `y <- detrend(x, tt = "linear")`, then `x` equals `y` + linear trend or
#' linear trend = `x - y`.
#' @export
#'
#' @examples
#' x <- 1:11
#' x - detrend(x, tt = "constant")
#' detrend(x, tt = "linear")
detrendcl <- function(x, tt) {
  if (!is.numeric(x) && !is.complex(x))
    stop("'x' must be a numeric or complex vector or matrix.")

  if (is.vector(x))
    x <- as.matrix(x)
  n <- nrow(x)

  if (tt == "constant") {
    y <- x - mean(x)
  }
  else if (tt == "linear") {
    a <- cbind(matrix(1, nrow = n, ncol = 1), matrix(1:n, nrow = n, ncol = 1) / n)
    # Get residuals from linear regression
    y <- x - a %*% solve(t(a) %*% a, t(a) %*% x)
  }
  else {
    stop("Trend type 'tt' must be 'constant' or 'linear'.")
  }
  return(as.vector(y))
}
