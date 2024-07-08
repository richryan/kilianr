## code to prepare `DATASET` dataset goes here

library(here)
library(readr)
library(dplyr)
library(lubridate)

fedfunds <- read_tsv(
  here("inst", "extdata", "kilian-lutkepohl-ch02_fedfunds.txt"),
  col_names = c("year", "month", "fedfunds"),
  col_types = "iid"
) |>
  mutate(date = ymd(paste(year, month, "01", sep = "-")),
         qtr = quarter(date, with_year = TRUE)) %>%
  group_by(qtr) %>%
  summarise(irate = mean(fedfunds))

gnp <- read_tsv(
  here("inst", "extdata", "kilian-lutkepohl-ch02_realgnp.txt"),
  col_names = c("year", "quarter", "rgnp"),
  col_types = "iid"
) |>
  mutate(date = yq(paste0(year, ":Q", quarter)),
         qtr = quarter(date, with_year = TRUE),
         drgnp = 100*(log(rgnp) - lag(log(rgnp), n = 1)))

deflator <- read_tsv(
  here("inst", "extdata", "kilian-lutkepohl-ch02_gnpdeflator.txt"),
  col_names = c("year", "quarter", "gnp_deflator"),
  col_types = "iid"
) |>
  mutate(date = yq(paste(year, ":Q", quarter)),
         qtr = quarter(date, with_year = TRUE),
         infl = 100*(log(gnp_deflator) - lag(log(gnp_deflator), n = 1)))

dat <- deflator %>%
  left_join(gnp, by = c("qtr", "year", "quarter", "date")) %>%
  left_join(fedfunds, by = c("qtr"))

# Keep data from 1954q4-2007q4
kilianLutkepohlCh02Macro <- dat %>%
  filter(qtr >= "1954.4") |>
  select(year, quarter, qtr, date, gnp_deflator, infl, rgnp, drgnp, irate)

usethis::use_data(kilianLutkepohlCh02Macro, overwrite = TRUE)
