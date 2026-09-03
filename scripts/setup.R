# 00_setup.R - packages, shared paths, shared settings, shared helpers

pkgs <- c("forecast", "tseries", "dplyr", "ggplot2", "lubridate", "gridExtra")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)
library(tseries)
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)

# create output folders
dir.create("data", showWarnings = FALSE)
dir.create("output/plots/eda",        showWarnings = FALSE, recursive = TRUE)
dir.create("output/plots/models",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/plots/comparison", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

# path helpers
fig <- function(category, name) file.path("output/plots", category, name)
tbl <- function(name) file.path("output/tables", name)

# y-axis labels in millions
scale_y_millions <- function() {
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
}

# ---------------------------------------------------------------------
# shared settings
# ---------------------------------------------------------------------

# holdout size: one full seasonal cycle
H <- 12

# max lag for Ljung-Box / ACF checks
LAG_MAX <- 16

# load MCO-resolved series
load_series <- function() readRDS("data/ampang_monthly_full_resolved.rds")

# train/test split
split_series <- function(y, h = H) {
  list(train = head(y, length(y) - h), test = tail(y, h))
}

# ---------------------------------------------------------------------
# shared diagnostic helpers
# ---------------------------------------------------------------------

# count residual ACF lags outside the 95% bound
acf_out_of_bounds <- function(resid, lag.max = LAG_MAX) {
  r <- na.omit(resid)
  n <- length(r)
  ci <- 1.96 / sqrt(n)
  a <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# degrees of freedom used by the fitted ARMA orders
model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
}

# Ljung-Box test
lb_test <- function(fit, lag) {
  Box.test(residuals(fit), lag = lag, fitdf = model_fitdf(fit),
           type = "Ljung-Box")
}

# train/test overfitting check
gap_check <- function(acc_matrix) {
  mape_tr <- acc_matrix["Training set", "MAPE"]
  mape_te <- acc_matrix["Test set", "MAPE"]
  mase_tr <- acc_matrix["Training set", "MASE"]
  mase_te <- acc_matrix["Test set", "MASE"]
  rmse_tr <- acc_matrix["Training set", "RMSE"]
  rmse_te <- acc_matrix["Test set", "RMSE"]

  mape_gap   <- abs(mape_te - mape_tr) / mape_te
  mase_gap   <- abs(mase_te - mase_tr) / mase_te
  rmse_ratio <- rmse_te / rmse_tr

  list(mape_gap   = mape_gap,
       mase_gap   = mase_gap,
       rmse_ratio = rmse_ratio,
       rmse_train = rmse_tr,
       mase_train = mase_tr,
       direction  = if (mase_te > mase_tr) "test worse" else "train worse",
       within_10pct = (mase_gap <= 0.10),
       within_1_3x  = (rmse_ratio <= 1.3))
}

# one-line summary row per model
model_summary <- function(name, fit, fc, test) {
  a  <- accuracy(fc, test)
  g  <- gap_check(a)
  r  <- residuals(fit)
  lb12  <- lb_test(fit, 12)
  lbmax <- lb_test(fit, LAG_MAX)

  data.frame(
    model        = name,
    MAPE_test    = round(a["Test set", "MAPE"], 3),
    RMSE_test    = round(a["Test set", "RMSE"], 0),
    MAE_test     = round(a["Test set", "MAE"], 0),
    MASE_test    = round(a["Test set", "MASE"], 3),
    RMSE_train   = round(g$rmse_train, 0),
    MASE_train   = round(g$mase_train, 3),
    fitdf        = model_fitdf(fit),
    lb_pvalue_12 = round(lb12$p.value, 4),
    lb_pvalue_16 = round(lbmax$p.value, 4),
    lb_pass      = (lb12$p.value > 0.05) & (lbmax$p.value > 0.05),
    n_lags_out_12 = acf_out_of_bounds(r, lag.max = 12),
    n_lags_out_16 = acf_out_of_bounds(r, lag.max = LAG_MAX),
    mase_gap_pct = round(g$mase_gap * 100, 1),
    within_10pct = g$within_10pct,
    rmse_ratio   = round(g$rmse_ratio, 3),
    within_1_3x  = g$within_1_3x,
    direction    = g$direction,
    mape_gap_pct = round(g$mape_gap * 100, 1),
    stringsAsFactors = FALSE
  )
}