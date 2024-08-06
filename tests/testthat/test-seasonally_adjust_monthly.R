test_that("seasonal adjustment works", {
  # Use example from seasonal package
  m <- seas(AirPassengers)
  out_seasonal <- as.vector(t(final(m)))

  df <- tibble::tibble(
    y_nsa = as.vector(t(AirPassengers)),
    # see start(AirPassengers) and end(AirPassengers)
    date = seq(lubridate::ymd('1949-01-01'), lubridate::ymd('1960-12-01'), by='months'),
    y = seasonally_adjust_monthly(y_nsa, date)
  )

  expect_equal(df$y, out_seasonal)
})


test_that("seasonal adjustment works with missing values padded at start", {
  # Use example from seasonal package
  m <- seas(AirPassengers)
  out_seasonal <- as.vector(t(final(m)))

  df <- tibble::tibble(
    y_nsa = as.vector(t(AirPassengers)),
    # see start(AirPassengers) and end(AirPassengers)
    date = seq(lubridate::ymd('1949-01-01'), lubridate::ymd('1960-12-01'), by='months'),
  )

  df_na <- tibble::tibble(
    date = seq(lubridate::ymd('1948-06-01'), lubridate::ymd('1948-12-01'), by='months'),
  )

  df <- dplyr::bind_rows(df, df_na) |>
    dplyr::arrange(date)

  df <- df |>
    dplyr::mutate(y = seasonally_adjust_monthly(y_nsa, date)) |>
    tidyr::drop_na(y)

  expect_equal(df$y, out_seasonal)
})



