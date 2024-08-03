#' Remove the mean (constant) or linear trend (linear) from a vector.
#'
#' Missing values are ignored. For example, suppose the change in production is
#' the variable of interest. Generating the variable generates a missing value
#' in the first position. A linear trend using `detrendcl()` can still be used.
#' But missing values in the middle of the series will also be ignored, which is
#' probably not the desired outcome. A message is thrown when a missing value is
#' encountered.
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
  if (!is.numeric(x) && !is.complex(x)) {
    stop("'x' must be a numeric or complex vector or matrix.")
  }
  if (is.vector(x)) {
    x <- as.matrix(x)
  }

  N <- length(x)

  x_indxnon <- which(!is.na(x))

  xx <- x[x_indxnon, , drop = FALSE]
  n <- nrow(xx)

  if (n != N) {
    message("   *** Missing values detected.")
  }

  if (tt == "constant") {
    y <- xx - mean(xx)

    x[x_indxnon] <- y
    stopifnot(near(mean(x, na.rm = TRUE), 0))
    # To remove matrix-name convention: myvar[,1]
    return(as.vector(x))
  }
  else if (tt == "linear") {
    a <- cbind(matrix(1, nrow = n, ncol = 1), matrix(1:n,
                                                     nrow = n, ncol = 1)/n)
    y <- xx - a %*% solve(t(a) %*% a, t(a) %*% xx)

    x[x_indxnon] <- y
    return(as.vector(x))
  }
  else {
    stop("Trend type 'tt' must be 'constant' or 'linear'.")
  }
  return(as.vector(y))
}


# mydetrendcl <- function (x, tt) {
#   if (!is.numeric(x) && !is.complex(x)) {
#     stop("'x' must be a numeric or complex vector or matrix.")
#   }
#   if (is.vector(x)) {
#     x <- as.matrix(x)
#   }
#
#   x_indxnon <- which(!is.na(x))
#
#   xx <- x[x_indxnon, , drop = FALSE]
#   n <- nrow(xx)
#   if (tt == "constant") {
#     y <- xx - mean(xx)
#
#     x[x_indxnon] <- y
#     stopifnot(near(mean(x, na.rm = TRUE), 0))
#     return(x)
#   }
#   else if (tt == "linear") {
#     a <- cbind(matrix(1, nrow = n, ncol = 1), matrix(1:n,
#                                                      nrow = n, ncol = 1)/n)
#     y <- xx - a %*% solve(t(a) %*% a, t(a) %*% xx)
#
#     x[x_indxnon] <- y
#     return(x)
#   }
#   else {
#     stop("Trend type 'tt' must be 'constant' or 'linear'.")
#   }
#   return(as.vector(y))
# }
#
# nx <- 33
# x <- 1:nx
# y <- 3 + 2 * x + rnorm(nx)
#
# z <- mydetrendcl(y, tt = "linear")
#
# plot(x, y)
# points(x, y - z, col = "blue", type = "l")
#
#
#

