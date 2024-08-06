test_that("seasonal adjustment works", {

  mytime <- seq(lubridate::ymd('1976-01-01'), lubridate::ymd('2020-12-01'), by='months')
  N <- length(mytime)
  y <- 0 + cos(seq(1:N) / 3)

  df <- tibble::tibble(
    date = mytime,
    y_nsa = y,
    y_sa = seasonally_adjust_monthly(y_nsa, date)
  )

  # ggplot2::ggplot(data = df) +
  #   ggplot2::geom_line(mapping = ggplot2::aes(x = date, y = y))

})


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
