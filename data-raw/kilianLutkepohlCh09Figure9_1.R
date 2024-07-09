## code to prepare `kilianLutkepohlCh09Figure9_1` dataset goes here
# Code that goes with Kilian and Lutkepohl's textbook,
# is provided at Kilian's website: https://sites.google.com/site/lkilian2019/textbook/code
# These notes were included with the data:
#
#    Quarterly data from FRED
#    GDP: Chain-type Price Index: Index 2009
#    Real Gross Domestic Product: Billions of Chained 2009 Dollars
#    1973.I-2013.II
#
#    Monthly WTI spot price from Economagic, 1973.1-2007.12
#
# I incorporated the data in the repository. This code converts the monthly data
# to a quarterly frequency by taking data from the last month of each quarter.
# This is what the original MATLAB code does. The oil dataset was lightly edited
# by me to include a tab in the first row.

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

dat_gdpdeflator <- read_tsv(
  here(
    "inst",
    "extdata",
    "kilian-lutkepohl-ch09_gdpdeflator2.txt"
  ),
  col_names = c("year", "quarter", "gdpdeflator"),
  col_types = "iid"
) |>
  mutate(lgdpdeflator = log(gdpdeflator),
         infl = 100 * (lgdpdeflator - lag(lgdpdeflator)))

dat_realgdp <- read_tsv(
  here("inst",
       "extdata",
       "kilian-lutkepohl-ch09_realgdp2.txt"),
  col_names = c("year", "quarter", "rgdp"),
  col_types = "iid"
) %>%
  mutate(lrgdp = log(rgdp),
         drgdp = 100 * (lrgdp - lag(lrgdp)))

# The original dataset did not have a tab in the first row,
# so this uses a lightly edited version of the file
dat_poil <- read_tsv(
  here("inst",
       "extdata",
       "kilian-lutkepohl-ch09_poil_edited.txt"),
  col_names = c("year", "month", "poil"),
  col_types = "iid"
) %>%
  mutate(date = ymd(paste(year, month, "01", sep = "-")),
         quarter = quarter(date, with_year = FALSE),
         lpoil = log(poil),
         dpoil = 100 * (lpoil - lag(lpoil, n = 3))) %>%
  # How the Kilian and Lutkepohl code goes from monthly to quarterly
  filter(month %in% c(3, 6, 9, 12))

dat <- dat_poil %>%
  left_join(dat_realgdp, by = c("year", "quarter")) %>%
  left_join(dat_gdpdeflator, by = c("year", "quarter")) %>%
  tidyr::drop_na(dpoil) %>%
  mutate(drpoil = dpoil - infl)

kilianLutkepohlCh09Figure9_1 <- dat |>
  select(year, quarter, date, drpoil, infl, drgdp)

# filter(year > 1986)
usethis::use_data(kilianLutkepohlCh09Figure9_1, overwrite = TRUE)
