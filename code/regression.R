# ==============================================================================
# Empirical Econometric Analysis: Climate Anomalies & Agricultural Returns
# Estimates OLS (Conditional Mean) and Quantile Regressions across percentiles:
# tau in {0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95}
# ==============================================================================

libraries <- c(
  "data.table", "zoo", "countrycode", "lubridate", "xts",
  "tidyverse", "plm", "stringr", "broom", "stargazer", "quantreg"
)

lapply(libraries, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
})

set.seed(123)

dataset_path <- file.path("data", "enso_shock_dataset.csv")
if (!file.exists(dataset_path)) {
  stop("Input dataset not found at: ", dataset_path, ". Please review data/README.md for reproduction instructions.")
}

enso_shock_dataset <- fread(dataset_path)
dir.create("results", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Baseline Ordinary Least Squares (OLS) - Conditional Mean
# ------------------------------------------------------------------------------
base_OLS_corn <- lm(d.corn ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset)
base_OLS_soybean <- lm(d.soybean ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset)
base_OLS_wheat <- lm(d.wheat ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset)

stargazer(
  list(base_OLS_corn, base_OLS_soybean, base_OLS_wheat),
  type = "html",
  out = "results/OLS_base.html",
  style = "aer",
  column.labels = c("d.corn", "d.soybean", "d.wheat")
)

# ------------------------------------------------------------------------------
# 2. Quantile Regressions across Tau Percentiles: tau in {0.05, 0.10, ..., 0.95}
# ------------------------------------------------------------------------------
taus <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

# Corn
qr_corn <- lapply(taus, function(tau) {
  rq(d.corn ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset, tau = tau)
})
stargazer(
  qr_corn,
  column.labels = as.character(taus),
  type = "html",
  out = "results/corn_qr_base.html",
  rq.se = "iid",
  style = "aer",
  initial.zero = FALSE,
  single.row = TRUE
)

# Soybean
qr_soybean <- lapply(taus, function(tau) {
  rq(d.soybean ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset, tau = tau)
})
stargazer(
  qr_soybean,
  column.labels = as.character(taus),
  type = "html",
  out = "results/soybean_qr_base.html",
  rq.se = "iid",
  style = "aer",
  initial.zero = FALSE,
  single.row = TRUE
)

# Wheat
qr_wheat <- lapply(taus, function(tau) {
  rq(d.wheat ~ `ANOM1+2` + d.sp500 + ads + d.gepui + d.vix + d.oil + shadow_rate + d.twfx, data = enso_shock_dataset, tau = tau)
})
stargazer(
  qr_wheat,
  column.labels = as.character(taus),
  type = "html",
  out = "results/wheat_qr_base.html",
  rq.se = "iid",
  style = "aer",
  initial.zero = FALSE,
  single.row = TRUE
)

message("Regressions executed successfully. Results written to results/ directory.")
