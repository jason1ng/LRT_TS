# 03_eda.R - EDA + stationarity/seasonality checks on the MCO-resolved,
# full-period (91-month) series. Each check here is what justifies later
# choosing SARIMA and its p,d,q / P,D,Q settings - not just plots for
# their own sake. Uses base R graphics for full layout/axis control.

source("scripts/setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

# helper: y-axis in millions instead of scientific notation
axis_millions <- function() {
  at <- pretty(par("usr")[3:4])
  axis(2, at = at, labels = paste0(round(at / 1e6, 1), "M"))
}

png(fig("eda", "acf_pacf_raw.png"), width = 900, height = 400, res = 130)
par(mfrow = c(1, 2))
Acf(ampang_ts,  lag.max = 24, main = "ACF (original series)")
Pacf(ampang_ts, lag.max = 24, main = "PACF (original series)")
dev.off()

# --- Visual EDA: combined into ONE figure (2x2 grid) ---
# Panel 1: raw series - shows trend + visible disruption/recovery shape
# Panel 2: seasonal plot - checks if a repeating annual pattern exists
# Panel 3/4: ACF/PACF (on differenced series) - checks for a lag-12 spike,
# which is the direct evidence for using SEASONAL (SARIMA) terms

png(fig("eda", "eda_summary.png"), width = 1000, height = 700, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# panel 1: raw series
ts.plot(ampang_ts, main = "Monthly Ridership (2019-2026, MCO-resolved)",
        xlab = "Year", ylab = "", yaxt = "n")
axis_millions()

# panel 2: seasonal plot (base graphics version of ggseasonplot)
seasonplot(ampang_ts, main = "Seasonal Plot", ylab = "", yaxt = "n",
           col = rainbow(8), year.labels = FALSE)
axis_millions()

# panel 3: ACF on first difference
Acf(diff(ampang_ts), lag.max = 24, main = "ACF (first difference)")

# panel 4: PACF on first difference
Pacf(diff(ampang_ts), lag.max = 24, main = "PACF (first difference)")

dev.off()

# STL decomposition - splits series into trend/seasonal/remainder, used
# to check seasonal shape/consistency and compute seasonal strength below
decomp <- stl(ampang_ts, s.window = "periodic", robust = TRUE)
png(fig("eda", "stl_decomposition.png"), width = 800, height = 600, res = 130)
plot(decomp, main = "STL Decomposition")
dev.off()

# seasonal plot, standalone full-size version
png(fig("eda", "seasonal_plot.png"), width = 800, height = 600, res = 130)
par(mar = c(4, 4, 3, 1))
seasonplot(ampang_ts, main = "Seasonal Plot", ylab = "", yaxt = "n",
           col = rainbow(8), year.labels = TRUE)
axis_millions()
dev.off()

# lag plot - series against itself at lags 1-12. A tight positive
# pattern specifically at lag 12 is the visual confirmation of annual
# seasonality, matching the ACF/PACF spike above
png(fig("eda", "lag_plot.png"), width = 800, height = 800, res = 130)
par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
for (k in 1:12) {
  d <- ts.intersect(y = ampang_ts, y_lag = stats::lag(ampang_ts, -k))
  plot(d[, "y_lag"], d[, "y"], xlab = paste0("lag ", k), ylab = "",
       main = paste0("lag ", k), pch = 20, col = "steelblue")
}
dev.off()

# seasonal strength (Hyndman's measure) - quantifies how much of the
# variation is explained by the seasonal pattern vs noise; the main
# numeric justification for choosing SARIMA over plain ARIMA
seasonal_strength <- max(0, 1 - var(decomp$time.series[, "remainder"]) /
                            var(decomp$time.series[, "seasonal"] + decomp$time.series[, "remainder"]))
cat("Seasonal strength:", round(seasonal_strength, 3), "\n")

# --- stationarity tests: determine differencing order (d) ---
cat("\n--- Level series ---\n")
print(adf.test(ampang_ts))   # want p < 0.05 for stationary
print(kpss.test(ampang_ts))  # want p > 0.05 for stationary

cat("\n--- First difference ---\n")
print(adf.test(diff(ampang_ts)))
print(kpss.test(diff(ampang_ts)))