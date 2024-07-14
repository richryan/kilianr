## code to prepare `kilianLutkepohlCh12Figure12_5` dataset goes here

library(here)
library(readr)
library(dplyr)
library(lubridate)

dat <- read_table(here("inst",
                       "extdata",
                       "kilian-lutkepohl-ch12_figure12-5.txt"),
                  col_names = c("oil_supply", "agg_demand", "rpoil")) |>
  mutate(date = seq.Date(
    from = ymd("1973-02-01"),
    to = ymd("2007-12-01"),
    by = "month"
  ))

kilianLutkepohlCh12Figure12_5 <- dat |>
  select(date, everything())

usethis::use_data(kilianLutkepohlCh12Figure12_5, overwrite = TRUE)
