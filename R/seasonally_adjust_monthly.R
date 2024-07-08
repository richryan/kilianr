#' Seasonally adjust monthly data
#'
#' `seasonally_adjust_monthly()` uses the function `seas` from the `seasonal` package,
#' which "calls the automatic procedures of X-13ARIMA-SEATS to perform a seasonal adjustment that works well in most circumstances."
#'
#' @param x A vector.
#' @param date A vector of monthly dates.
#'
#' @returns A vector.
#' @importFrom seasonal seas
#' @export
#'
#' @examples
#' seasonally_adjust(data$some_series, data$date)
seasonally_adjust_monthly <- function(x, date) {
  junk <- final(seasonal::seas(ts(
    x,
    start = c(year(min(date)), month(min(date))),
    frequency = 12
  )))
  return(as.numeric(junk))
}
