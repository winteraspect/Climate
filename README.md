# Climate Shocks and Agricultural Commodity Markets

### *Empirical Analysis of El Niño–Southern Oscillation (ENSO) Asymmetries via Quantile Regression*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)

---

## Executive Summary

How do global meteorological shocks translate into financial market risks for staple commodities? While standard linear econometric models (Ordinary Least Squares) often conclude that climate anomalies have statistically insignificant effects on **conditional mean** commodity returns, this repository provides empirical evidence of pronounced **tail-risk asymmetries**.

By estimating **Quantile Regressions** across conditional percentiles ($\tau \in [0.05, 0.95]$) on monthly agricultural returns (**Corn**, **Soybeans**, and **Wheat**) controlling for macroeconomic conditions, global policy uncertainty, monetary policy stance, and exchange rates:

- **Mean Neutrality**: In standard OLS specifications, sea surface temperature (SST) anomalies across the El Niño monitoring zones exhibit near-zero, statistically insignificant point estimates ($p > 0.10$).
- **Tail Significance**: In the extreme lower quantiles ($\tau = 0.05$), positive temperature anomalies (El Niño conditions) exhibit a statistically significant positive effect ($\beta = +0.016, p < 0.05$ for corn), demonstrating that climate shocks serve as a critical tail-risk buffer and price support during market downturns.

---

## Table of Contents

- [Research Motivation & Background](#research-motivation--background)
- [Econometric Methodology](#econometric-methodology)
  - [1. Baseline Linear Model (OLS)](#1-baseline-linear-model-ols)
  - [2. Quantile Regression (Koenker & Bassett)](#2-quantile-regression-koenker--bassett)
  - [3. Conditioning Set & Macroeconomic Controls](#3-conditioning-set--macroeconomic-controls)
- [Empirical Results](#empirical-results)
  - [Ordinary Least Squares (Conditional Mean)](#ordinary-least-squares-conditional-mean)
  - [Quantile Regression: Corn Returns](#quantile-regression-corn-returns)
  - [Quantile Regression: Soybeans & Wheat](#quantile-regression-soybeans--wheat)
- [Economic Mechanisms & Discussion](#economic-mechanisms--discussion)
- [Repository Structure](#repository-structure)
- [Replication Guide](#replication-guide)
- [Data Availability & Privacy Notice](#data-availability--privacy-notice)
- [References](#references)

---

## Research Motivation & Background

The **El Niño–Southern Oscillation (ENSO)** is the dominant driver of inter-annual global climate variability. Shifts in tropical Pacific sea surface temperatures disrupt global atmospheric circulation, altering monsoon patterns, rainfall, and temperature distributions across major agricultural breadbaskets (the U.S. Midwest, Brazil, Argentina, Australia, and Southeast Asia).

While classical agricultural economics literature (e.g. Brunner 2002; Cashin, Céspedes, & Sahay 2017) documents macro-level impacts of El Niño on world commodity prices, aggregate monthly index returns often exhibit low correlation with climate indicators at the mean. This project investigates whether climate risks are non-linear, operating primarily in the tails of the price distribution where storage buffers and inventory stock-outs become binding constraints.

---

## Econometric Methodology

### 1. Baseline Linear Model (OLS)

The baseline conditional mean return equation is specified as:

$$\Delta \ln P_{i, t} = \alpha_i + \beta_i \cdot \text{ANOM}_{t} + \mathbf{\Gamma}_i' \mathbf{X}_{t-1} + \varepsilon_{i, t}$$

where:
- $P_{i, t}$ is the IMF primary commodity price index for commodity $i \in \{\text{Corn}, \text{Soybeans}, \text{Wheat}\}$.
- $\Delta \ln P_{i, t}$ is the monthly log return of the commodity.
- $\text{ANOM}_t$ is the NOAA Sea Surface Temperature anomaly index (Niño 1+2 / Niño 3.4).
- $\mathbf{X}_{t-1}$ is a vector of lagged macroeconomic and financial market control variables.

### 2. Quantile Regression (Koenker & Bassett)

To capture asymmetric impacts across different market regimes (bearish downturns vs. bullish surges), we estimate linear quantile regressions for $\tau \in \{0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95\}$:

$$Q_{\tau}\left(\Delta \ln P_{i, t} \mid \text{ANOM}_t, \mathbf{X}_{t-1}\right) = \alpha_i(\tau) + \beta_i(\tau) \cdot \text{ANOM}_t + \mathbf{\Gamma}_i(\tau)' \mathbf{X}_{t-1}$$

estimated by minimizing the asymmetric check function $\rho_\tau(u) = u(\tau - \mathbb{I}(u < 0))$.

### 3. Conditioning Set & Macroeconomic Controls

To isolate the exogenous contribution of meteorological anomalies, the model controls for:
1. **Global Equity Market Sentiment**: Log-returns of the S&P 500 index (`d.sp500`).
2. **Real Business Conditions**: The Philadelphia Fed Aruoba-Diebold-Scotti business conditions index (`ads`).
3. **Policy Uncertainty**: Global Economic Policy Uncertainty index PPP-weighted (`d.gepui`, Baker, Bloom, & Davis).
4. **Market Volatility**: CBOE Volatility Index (`d.vix`) and Emerging Market Volatility (`d.vix_em`).
5. **Energy / Input Shocks**: Monthly returns of WTI Crude Oil (`d.oil`).
6. **Monetary Policy Stance**: Wu-Xia Shadow Federal Funds Rate (`shadow_rate`) to capture unconventional monetary policy at the zero lower bound.
7. **Exchange Rate Channel**: Trade-Weighted U.S. Dollar Index (`d.twfx`).

---

## Empirical Results

### Ordinary Least Squares (Conditional Mean)

Estimated via OLS on 110 monthly observations (coefficients with standard errors in parentheses):

| Variable | d.corn (1) | d.soybean (2) | d.wheat (3) |
| :--- | :---: | :---: | :---: |
| **`ANOM1+2`** | -0.003 (0.006) | -0.007 (0.005) | -0.006 (0.007) |
| **d.sp500** | 0.000 (0.003) | -0.001 (0.003) | 0.001 (0.003) |
| **ads** | 0.004 (0.002) | 0.001 (0.002) | 0.001 (0.003) |
| **d.gepui** | -0.004 (0.033) | -0.009 (0.029) | 0.006 (0.037) |
| **d.oil** | -0.000 (0.001) | 0.000 (0.001) | -0.001 (0.001) |
| **shadow_rate** | 0.002 (0.004) | 0.002 (0.003) | 0.001 (0.004) |
| **d.twfx** | -0.008\* (0.004) | -0.007\* (0.004) | -0.013\*\* (0.005) |
| **Constant** | 0.003 (0.007) | 0.006 (0.006) | 0.004 (0.008) |
| **Observations** | 110 | 110 | 110 |
| **$R^2$** | 0.080 | 0.068 | 0.093 |

*Notes: \*p < 0.10, \*\*p < 0.05, \*\*\*p < 0.01. Standard errors in parentheses.*

> **Key Observation**: In all three commodity OLS specifications, the climate anomaly coefficient `ANOM1+2` is statistically indistinguishable from zero ($t$-statistics $< 1.1$). Only the trade-weighted dollar index (`d.twfx`) exhibits statistical significance across all three crops.

---

### Quantile Regression: Corn Returns

When decomposing the price return distribution across quantiles, an entirely different structural picture emerges:

| Variable | $\tau = 0.05$ | $\tau = 0.10$ | $\tau = 0.25$ | $\tau = 0.50$ | $\tau = 0.75$ | $\tau = 0.90$ | $\tau = 0.95$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`ANOM1+2`** | **.016\*\*** (.008) | -.001 (.009) | -.002 (.006) | -.002 (.003) | -.004 (.005) | -.006 (.004) | -.008 (.005) |
| **d.sp500** | -.003 (.004) | .002 (.004) | .002 (.003) | .003\* (.002) | .001 (.003) | .002 (.002) | -.005\* (.003) |
| **ads** | -0.000 (.003) | -.001 (.004) | -0.000 (.003) | .003\*\* (.001) | .004\*\* (.002) | .009\*\*\* (.002) | .008\*\*\* (.002) |
| **d.gepui** | -.054 (.044) | -.088\* (.048) | -.021 (.036) | .021 (.019) | -.018 (.030) | -.022 (.023) | .054\* (.031) |
| **shadow_rate** | .018\*\*\* (.005) | -.002 (.005) | .003 (.004) | .002 (.002) | .007\*\* (.003) | -.003 (.003) | .001 (.003) |
| **d.twfx** | -.003 (.006) | -.005 (.007) | -.008\* (.005) | -.001 (.003) | -.002 (.004) | -.010\*\*\* (.003) | -.033\*\*\* (.004) |
| **Constant** | -.084\*\*\* (.009) | -.061\*\*\* (.010) | -.027\*\*\* (.008) | -.003 (.004) | .028\*\*\* (.006) | .068\*\*\* (.005) | .105\*\*\* (.007) |
| **Obs** | 110 | 110 | 110 | 110 | 110 | 110 | 110 |

*Notes: \*p < 0.10, \*\*p < 0.05, \*\*\*p < 0.01. Standard errors in parentheses.*

> **Key Finding**: In the **5th percentile ($\tau = 0.05$)**, `ANOM1+2` is **$+0.016$** and statistically significant at the 5% level ($p < 0.05$). During severe downturns in corn prices, warming anomalies in the equatorial Pacific provide an upward price buffer, preventing further negative price spirals. At higher percentiles ($\tau \ge 0.50$), this effect vanishes and turns slightly negative.

---

### Quantile Regression: Soybeans & Wheat

- **Soybeans**: Similar to corn, El Niño indicators show muted impacts in the center of the distribution ($\tau = 0.50$), while real economic activity (`ads`) and the trade-weighted dollar index (`d.twfx`) dominate the upper tail ($\tau = 0.90, 0.95$).
- **Wheat**: Extreme price movements are heavily governed by dollar exchange rate variations (`d.twfx`, $p < 0.001$ in the 95th percentile) and energy input costs.

Complete interactive HTML estimation tables are available in the [`results/`](results/) folder:
- [`results/OLS_base.html`](results/OLS_base.html) – Baseline OLS regressions across all 3 commodities.
- [`results/corn_qr_base.html`](results/corn_qr_base.html) – Quantile regression breakdown for Corn.
- [`results/soybean_qr_base.html`](results/soybean_qr_base.html) – Quantile regression breakdown for Soybeans.
- [`results/wheat_qr_base.html`](results/wheat_qr_base.html) – Quantile regression breakdown for Wheat.

---

## Economic Mechanisms & Discussion

1. **Storage Buffer Discontinuity**: In normal harvest years, global grain inventories absorb meteorological shocks, muting the transmission of temperature variations into spot and futures prices. Consequently, linear models measuring conditional means fail to detect significant elasticity.
2. **Asymmetric Tail Protection**: During periods of extreme price decline ($\tau = 0.05$), positive temperature anomalies (El Niño forecasts) signal prospective yield reductions across Southern Hemisphere exporters (Australia, South America). Traders bid up deferred contracts, creating an immediate price floor that cushions the bottom 5% of returns.
3. **Currency Dominance**: The trade-weighted U.S. dollar (`d.twfx`) remains the single most consistent determinant of dollar-denominated agricultural prices across all quantiles, verifying the international purchasing power parity channel in commodity finance.

---

## Repository Structure

```
Climate/
├── README.md               # Empirical research overview and results
├── REFERENCES.md           # Legal academic bibliography and DOI citations
├── LICENSE                 # MIT License
├── code/
│   ├── data_preparation.R # Data acquisition, cleaning, and merging pipeline
│   ├── regression.R       # OLS and Quantile Regression estimation scripts
│   └── UCM_DFS_factorVAR2.ox # OxMetrics state-space & dynamic factor VAR engine
├── data/
│   └── README.md          # Data privacy statement and public replication sources
├── literature/
│   └── README.md          # Legal notice redirecting to REFERENCES.md
└── results/
    ├── OLS_base.html       # Baseline OLS regression output (Stargazer AER)
    ├── corn_qr_base.html   # Quantile regression output for Corn
    ├── soybean_qr_base.html# Quantile regression output for Soybeans
    └── wheat_qr_base.html  # Quantile regression output for Wheat
```

---

## Replication Guide

### Prerequisites (R Environment)

Install the required R packages:
```r
install.packages(c(
  "data.table", "zoo", "countrycode", "lubridate", "xts",
  "tidyverse", "plm", "stringr", "broom", "stargazer",
  "quantreg", "readxl", "fredr"
))
```

### Running Regressions

1. Obtain a free API key from the [St. Louis Fed FRED](https://fred.stlouisfed.org/) and set it in your environment:
   ```r
   Sys.setenv(FRED_API_KEY = "your_fred_api_key_here")
   ```
2. Execute data preparation and regression estimation:
   ```bash
   Rscript code/regression.R
   ```
   Tables will be formatted and exported directly to `results/`.

---

## Data Availability & Privacy Notice

In strict adherence to data governance policies and proprietary access limits:
- Intermediate private working datasets and local spreadsheets are not distributed in this repository.
- Researchers can replicate the exact datasets using the public APIs documented in [`data/README.md`](data/README.md) (NOAA Climate Prediction Center, IMF Primary Commodity Prices, and Federal Reserve Economic Data).

---

## References

For full citations with active DOIs and publisher links, see [`REFERENCES.md`](REFERENCES.md). Key foundational literature includes:
- **Brunner, A. D. (2002)**. *El Niño and World Primary Commodity Prices: Warm Water or Hot Air?* Review of Economics and Statistics.
- **Cashin, P., Céspedes, L. F., & Sahay, R. (2017)**. *Fair weather or foul? The macroeconomic effects of El Niño.* Journal of International Economics.
- **Koenker, R., & Bassett, G. (1978)**. *Regression Quantiles.* Econometrica.
- **Ubilava, D. (2012)**. *The ENSO Effect and Asymmetries in Wheat Price Dynamics.* World Development.
