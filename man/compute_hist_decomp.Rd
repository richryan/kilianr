% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_hist_decomp.R
\name{compute_hist_decomp}
\alias{compute_hist_decomp}
\title{Compute historical decomposition}
\usage{
compute_hist_decomp(olsvarc_obj, date, var_order)
}
\arguments{
\item{olsvarc_obj}{Object created by \code{olsvarc()}}

\item{date}{The \code{date} from the original data (trimmed to T - p)}

\item{var_order}{Order of variables for structural VAR model}
}
\value{
A tibble dataframe that shows how
}
\description{
Compute historical decomposition
}
\examples{
y <- kilianLutkepohlCh12Figure12_5 |>
  dplyr::rename(oilsupply = oil_supply,
                aggdemand = agg_demand) |>
  dplyr::select(-date)
sol <- olsvarc(y, p = 24)
hd <- compute_hist_decomp(sol, date = kilianLutkepohlCh12Figure12_5$date, var_order = c("oilsupply", "aggdemand", "rpoil"))
plot(hd$date, hd$hd_series_shock_oilsupply_oilsupply, type = "l")
plot(hd$date, hd$hd_series_shock_aggdemand_oilsupply + hd$hd_series_shock_aggdemand_aggdemand + hd$hd_series_shock_aggdemand_rpoil, type = "l", col = "red")
points(kilianLutkepohlCh12Figure12_5$date, kilianLutkepohlCh12Figure12_5$agg_demand, type = "l")

# Another example ---------------------------------------------------------
y <- kilianLutkepohlCh02Macro |> select(drgnp, irate, infl)
m <- olsvarc(y, p = 4)
hd <- compute_hist_decomp(m, date = kilianLutkepohlCh02Macro$date, var_order = c("drgnp", "irate", "infl"))
}
