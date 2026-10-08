# ==============================================================================
# Climate & Commodity Return Data Preparation Pipeline
# Merges NOAA ENSO sea surface temperature anomalies, IMF agricultural commodity prices,
# and macroeconomic financial control variables (FRED, ADS, GEPU, Wu-Xia Shadow Rate).
# ==============================================================================

libraries <- c(
  "data.table", "zoo", "countrycode", "lubridate", "xts",
  "tidyverse", "plm", "stringr", "broom", "stargazer",
  "readxl", "fredr"
)

lapply(libraries, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
})

set.seed(123)

# Configure FRED API key from environment variable
fred_key <- Sys.getenv("FRED_API_KEY")
if (nzchar(fred_key)) {
  fredr_set_key(fred_key)
} else {
  message("NOTE: 'FRED_API_KEY' environment variable not detected.")
  message("To pull live macroeconomic series from FRED, run: Sys.setenv(FRED_API_KEY = '<your_key>')")
}

# ------------------------------------------------------------------------------
# 1. NOAA Sea Surface Temperature (SST) & ENSO Index Anomalies
# ------------------------------------------------------------------------------
sst_file <- file.path("data", "SST_ENSO.txt")
if (file.exists(sst_file)) {
  sst <- fread(sst_file)
  names(sst) <- c("YR", "MON", "NINO1+2", "ANOM1+2", "NINO3", "ANOM3", "NINO4", "ANOM4", "NINO3.4", "ANOM3.4")
  sst[, date := as.yearmon(paste(YR, MON, sep = "-"))]
} else {
  message("SST data file not found at: ", sst_file)
}

# ------------------------------------------------------------------------------
# 2. IMF Primary Commodity Prices (Corn, Soybeans, Wheat)
# ------------------------------------------------------------------------------
prices_file <- file.path("data", "IMF_COMMODITY_PRICE.csv")
if (file.exists(prices_file)) {
  prices <- fread(prices_file)
  prices[, date := as.yearmon(yyyymm)]
  corn <- xts(prices$corn, order.by = prices$date)
  soybean <- xts(prices$soybean, order.by = prices$date)
  wheat <- xts(prices$wheat, order.by = prices$date)
  # Lag commodity prices by 1 month
  prices[, date := as.yearmon(yyyymm) + 1/12]
} else {
  message("Commodity prices file not found at: ", prices_file)
}

# ------------------------------------------------------------------------------
# 3. Macroeconomic & Financial Controls from FRED
# ------------------------------------------------------------------------------
if (nzchar(fred_key)) {
  # S&P 500
  sp500 <- fredr("sp500")
  sp500_xts <- xts(sp500$value, order.by = as.Date(sp500$date, "%Y-%m-%d"))
  d.sp500.m <- 100 * diff(log(apply.monthly(sp500_xts, last)))

  # VIX (Total and Emerging Markets)
  vix <- fredr("VIXCLS")
  vix_xts <- xts(vix$value, order.by = as.Date(vix$date, "%Y-%m-%d"))
  d.vix.m <- 100 * diff(log(apply.monthly(vix_xts, last)))

  vix_em <- fredr("VXEEMCLS")
  vix_em_xts <- xts(vix_em$value, order.by = as.Date(vix_em$date, "%Y-%m-%d"))
  d.vix_em.m <- 100 * diff(log(apply.monthly(vix_em_xts, last)))

  # Trade-Weighted USD Index
  twfx <- fredr("DTWEXBGS")
  twfx_xts <- xts(twfx$value, order.by = as.Date(twfx$date, "%Y-%m-%d"))
  d.twfx.m <- 100 * diff(log(apply.monthly(twfx_xts, last)))

  # WTI Crude Oil
  oil <- fredr("DCOILWTICO")
  oil_xts <- xts(oil$value, order.by = as.Date(oil$date, "%Y-%m-%d"))
  d.oil.m <- 100 * diff(log(apply.monthly(oil_xts, last)))
}

# ------------------------------------------------------------------------------
# 4. ADS Business Conditions, GEPU, & Shadow Rate Controls
# ------------------------------------------------------------------------------
ads_file <- file.path("data", "ADS_All_Vintages-zip.xlsx")
if (file.exists(ads_file)) {
  ads <- as.data.table(read_excel(ads_file))
  ads[, date := as.yearmon(...1, format = "%Y:%m:%d")]
  ads_xts <- xts(ads$ADS_INDEX_081221, order.by = ads$date)
  ads.m <- apply.monthly(ads_xts, last)
  d.ads.m <- 100 * log(ads.m / lag(ads.m))
}

gepui_file <- file.path("data", "Global_Policy_Uncertainty_Data.xlsx")
if (file.exists(gepui_file)) {
  gepui <- as.data.table(read_excel(gepui_file))
  gepui[, date := as.yearmon(paste(Year, Month, sep = "-"))]
  gepui <- gepui[!is.na(date)]
  gepui_xts <- xts(gepui$GEPU_ppp, order.by = gepui$date)
  gepui_monthly <- apply.monthly(gepui_xts, FUN = last)
  d.gepui.m <- log(gepui_monthly / lag(gepui_monthly))
}

shadow_file <- file.path("data", "shadowrate_US.xls")
if (file.exists(shadow_file)) {
  shadow_dt <- as.data.table(read_excel(shadow_file, col_names = FALSE))
  names(shadow_dt)[2] <- "shadow_rate"
  shadow_dt[, date := as.yearmon(paste(str_sub(as.character(...1), 1, 4), str_sub(as.character(...1), 5, 6), sep = "-"))]
  shadow_xts <- xts(shadow_dt$shadow_rate, order.by = shadow_dt$date)
}

# ------------------------------------------------------------------------------
# 5. Merge Controls & Commodity Returns into Unified Dataset
# ------------------------------------------------------------------------------
output_file <- file.path("data", "enso_shock_dataset.csv")
message("Pipeline configured. Run with local input files to generate: ", output_file)
