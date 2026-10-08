# Data Availability & Sources

## Data Privacy & Governance Statement
In accordance with data governance policies and licensing constraints, raw intermediate datasets and proprietary extracts are not hosted directly in this public repository.

Researchers and practitioners can reproduce all datasets used in this project through the public APIs and open data repositories detailed below.

---

## Public Data Sources

### 1. El Niño-Southern Oscillation (ENSO) & Sea Surface Temperature (SST) Anomalies
- **Source**: National Oceanic and Atmospheric Administration (NOAA) – Climate Prediction Center (CPC)
- **Variables**: Niño 1+2, Niño 3, Niño 4, and Niño 3.4 monthly Sea Surface Temperature index and anomalies (base period 1991–2020).
- **Access**: [NOAA CPC Monthly Atmospheric & SST Indices](https://www.cpc.ncep.noaa.gov/data/indices/)

### 2. Primary Agricultural Commodity Prices
- **Source**: International Monetary Fund (IMF) Primary Commodity Price System (PCPS)
- **Variables**: Monthly price indices (USD/metric ton or index level) for:
  - Corn (Maize)
  - Soybeans
  - Wheat
- **Access**: [IMF Primary Commodity Prices Portal](https://www.imf.org/en/Research/commodity-prices)

### 3. Macroeconomic & Financial Controls
- **Source**: Federal Reserve Bank of St. Louis (FRED) API
- **Series IDs**:
  - `SP500`: S&P 500 Stock Price Index
  - `VIXCLS`: CBOE Volatility Index (VIX)
  - `VXEEMCLS`: CBOE Emerging Markets ETF Volatility Index
  - `DCOILWTICO`: Crude Oil Prices: West Texas Intermediate (WTI)
  - `DTWEXBGS`: Trade-Weighted U.S. Dollar Index (Broad Goods and Services)
- **Access**: [FRED Economic Data](https://fred.stlouisfed.org/) (requires free API key via `Sys.setenv(FRED_API_KEY = "...")`).

### 4. Real-Time Business Conditions
- **Source**: Federal Reserve Bank of Philadelphia (Aruoba, Diebold, & Scotti)
- **Series**: ADS Business Conditions Index (tracking real business conditions at high frequency).
- **Access**: [Philadelphia Fed ADS Index](https://www.philadelphiafed.org/surveys-data/real-time-data-research/ads)

### 5. Economic Policy Uncertainty
- **Source**: Economic Policy Uncertainty Index (Baker, Bloom, & Davis)
- **Series**: Global Economic Policy Uncertainty (GEPU) PPP-weighted index.
- **Access**: [PolicyUncertainty.com](https://www.policyuncertainty.com/)

### 6. Shadow Federal Funds Rate (Monetary Policy Stance)
- **Source**: Cynthia Wu & Fan Dora Xia (Wu-Xia Shadow Rate)
- **Access**: [Federal Reserve Bank of Atlanta Shadow Rate Data](https://www.atlantafed.org/cqer/research/shadow_rate)
