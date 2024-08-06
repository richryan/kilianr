#' Seasonally adjust monthly data
#'
#' `seasonally_adjust_monthly()` uses the function `seas` from the `seasonal` package,
#' which "calls the automatic procedures of X-13ARIMA-SEATS to perform a seasonal adjustment that works well in most circumstances."
#'
#' @param x A vector.
#' @param date A vector of monthly dates.
#'
#' @returns A vector.
#' @import seasonal
#' @import lubridate
#' @export
#'
#' @examples
#' seasonally_adjust(data$some_series, data$date)
seasonally_adjust_monthly <- function(x, date) {

  N <- length(x)

  x_indxnon <- which(!is.na(x))

  xx <- x[x_indxnon]
  n <- length(xx)

  if (n != N) {
    message("   *** Missing values detected.")
  }

  date_xx <- date[x_indxnon]

  junk <- seasonal::final(seasonal::seas(ts(
    xx,
    start = c(year(min(date_xx)), month(min(date_xx))),
    frequency = 12
  )))

  x[x_indxnon] <- as.numeric(junk)

  return(x)
}
