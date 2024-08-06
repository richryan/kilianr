# compute_hist_decomp <- function(variables) {
#
# }
#
# library(tidyr)
# library(dplyr)
y_global_oil <- kilianLutkepohlCh12Figure12_5 |>
  # mutate(check_oil_supply = detrendcl(oil_supply, tt = "constant")) |>
  mutate(across(c(oil_supply, agg_demand, rpoil), \(x) detrendcl(x, tt = "constant"))) |>
  rename(oilsupply = oil_supply,
         aggdemand = agg_demand) |>
  select(-date) |>
  as.matrix()

horizon <- 15

sol_global_oil <- olsvarc(y_global_oil, p = 24)
var_order_global_oil <- c("oilsupply", "aggdemand", "rpoil")
var_cumsum_global_oil <- c("response_shock_oilsupply_oilsupply",
                           "response_shock_oilsupply_aggdemand",
                           "response_shock_oilsupply_rpoil")
negative_shocks_global_oil <- c("oilsupply")

B0inv <- t(chol(sol_global_oil$SIGMAhat))

irf_global_oil <- irfvar(Ahat = sol_global_oil$Ahat,
                         B0inv = B0inv,
                         p = sol_global_oil$p, h = nrow(sol_global_oil$y) - sol_global_oil$p - 1,
                         negative_shocks = negative_shocks_global_oil,
                         var_cumsum = var_cumsum_global_oil,
                         var_order = var_order_global_oil)

IRF <- irf_global_oil$irfm

Ehat <- B0inv %*% sol_global_oil$Uhat

dim(IRF)
dim(Ehat)
dim(IRF) - dim(Ehat)

# # Keep the tidy version of responses
# dat_irf_global_oil <- irf_global_oil$irf_tidy
