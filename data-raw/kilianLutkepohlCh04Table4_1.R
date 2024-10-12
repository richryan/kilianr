## code to prepare `kilianLutkepohlCh04Table4_1` dataset goes here
# Code that goes with Kilian and Lutkepohl's textbook,
# is provided at Kilian's website: https://sites.google.com/site/lkilian2019/textbook/code
#
# The data from
#   Kilian, L. and Park, C. (2009), THE IMPACT OF OIL PRICE SHOCKS ON THE U.S.
#   STOCK MARKET. International Economic Review, 50: 1267-1287.
#   https://doi.org/10.1111/j.1468-2354.2009.00568.x
# and span 1973.2-2006.12, where
#
#   drpod: Global oil production growth
#   rea: Kilian (2009, AER) business cycle index of global real economic activity
#   rpoil: Real price of oil in percent deviations from mean
#   dd: CRSP real dividend growth in percent
#
# I incorporated the data in the repository.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)


dat_raw <- read_table(
  here(
    "inst",
    "extdata",
    "kilian-lutkepohl-ch04_table4-1.txt"
  ),
  col_names = c("dprod", "rea", "rpoil", "dd"),
  col_types = "dddd"
)

kilianLutkepohlCh04Table4_1 <- dat_raw |>
  mutate(date = seq(ymd("1973-02-01"), ymd("2006-12-01"), by = "months"))

usethis::use_data(kilianLutkepohlCh04Table4_1, overwrite = TRUE)
