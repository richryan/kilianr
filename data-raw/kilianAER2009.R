## code to prepare `kilianAER2009` dataset goes here
library(here)
library(readr)
library(dplyr)
library(lubridate)

fin <- here::here("inst", "extdata", "data.txt")
dat_raw <- readr::read_table(fin, col_names = FALSE)

date_start <- lubridate::ymd("1973-02-01")

kilian2009AER <- dat_raw |>
  dplyr::mutate(d_production = as.numeric(X1),
         real_economic_activity = as.numeric(X2),
         real_price_oil = as.numeric(X3),
         junk = row_number() - 1,
         date = date_start + months(junk)) |>
  select(d_production, real_economic_activity, real_price_oil, date)

usethis::use_data(kilian2009AER, overwrite = TRUE)
