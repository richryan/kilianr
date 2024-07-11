#' Data used in Lutz Kilian's "Not All Oil Price Shocks Are Alike: Disentangling Demand and Supply Shocks in the Crude Oil Market," American Economic Review, 99 (3): 1053--69.
#'
#' @format ## `kilian2009AER`
#' A data frame with 419 rows (data from 1973-02-01 to 2007-12-01) and 4 columns:
#' \describe{
#'   \item{d_production}{Change in oil production}
#'   \item{real_economic_activity}{Real economic activity}
#'   \item{real_price_oil}{Real price of oil}
#'   \item{date}{Date object (constructed with `lubridate::ymd()`)}
#' }
#' @source <https://www.aeaweb.org/articles?id=10.1257/aer.99.3.1053>
"kilian2009AER"

#' Data used in section 2.3, figure 2.1 of Kilian and Lutkepohl's book "Structural Vector Autoregressive Analysis."
#'
#' The data are US macro data on real GNP, the federal funds rate, and inflation.
#' The data are quarterly.
#'
#' @format ## `kilianLutkepohlCh02Macro`
#' A data frame with 213 rows (1954 quarter 4 to 2007 quarter 4) and 9 columns:
#' \describe{
#'   \item{year}{Year (integer)}
#'   \item{quarter}{Quarter (integer)}
#'   \item{qtr}{Quarter (Date object constructed with `lubridate::yq()`)}
#'   \item{date}{Date object}
#'   \item{gnp_deflator}{GNP deflator}
#'   \item{infl}{Inflation measured by 100 times the difference in the log of the GNP deflator}
#'   \item{rgnp}{Real economic activity}
#'   \item{drgnp}{Difference in log of real GNP times 100}
#'   \item{irate}{Quarterly average of monthly values of the federal funds rate}
#' }
#' @source Lutz Kilian's website: <https://sites.google.com/site/lkilian2019/textbook/code>
"kilianLutkepohlCh02Macro"

#' Data used in section 12.13, figure 12.5 of Kilian and Lutkepohl's book "Structural Vector Autoregressive Analysis."
#'
#' The data are monthly data on world crude-oil production, global real economic activity, and the real price of oil.
#'
#' @format ## `kilianLutkepohlCh12Figure12_5`
#' A data frame with 419 rows (data from 1973-02-01 to 2007-12-01) and 4 columns:
#' \describe{
#'   \item{date}{Date object, monthly (constructed with `lubridate::ymd()`)}
#'   \item{oil_supply}{Percent change in world crude-oil production}
#'   \item{agg_demand}{Index of real economic activity}
#'   \item{rpoil}{Real price of oil}
#' }
#' @source Lutz Kilian's website: <https://sites.google.com/site/lkilian2019/textbook/code>
"kilianLutkepohlCh12Figure12_5"

#' Data used in section 12.13, figure 12.5 of Kilian and Lutkepohl's book "Structural Vector Autoregressive Analysis."
#'
#' The data are monthly data on world crude-oil production, global real economic activity, and the real price of oil.
#'
#' @format ## `kilianLutkepohlCh09Figure9_1`
#' A data frame with 162 rows (data from 1973q1 to 2013q2) and 6 columns:
#' \describe{
#'   \item{year}{Year}
#'   \item{quarter}{Quarter}
#'   \item{date}{Date object, monthly (constructed with `lubridate::ymd()`)}
#'   \item{drpoil}{Percent change in the real WTI price of crude oil}
#'   \item{infl}{Inflation rate measured by the US GDP deflator}
#'   \item{drgdp}{Growth rate of US real GDP}
#' }
#' @source Lutz Kilian's website: <https://sites.google.com/site/lkilian2019/textbook/code>
"kilianLutkepohlCh12Figure12_5"
