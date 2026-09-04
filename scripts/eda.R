# Exploratory analysis and stationarity checks for the resolved series.

source("scripts/setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

# Label the y-axis in millions.
axis_millions <- function() {
  at <- pretty(par("usr")[3:4])
  axis(2, at = at, labels = paste0(round(at / 1e6, 1), "M"))
}

png(fig("eda", "acf_pacf_raw.png"), width = 900, height = 400, res = 130)
par(mfrow = c(1, 2))
Acf(ampang_ts,  lag.max = 24, main = "ACF (original series)")
Pacf(ampang_ts, lag.max = 24, main = "PACF (original series)")
dev.off()

# Summary figure: series, seasonal pattern, ACF, and PACF.

png(fig("eda", "eda_summary.png"), width = 1000, height = 700, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Original series.
ts.plot(ampang_ts, main = "Monthly Ridership (2019-2026, MCO-resolved)",
        xlab = "Year", ylab = "", yaxt = "n")
axis_millions()

# Seasonal pattern.
seasonplot(ampang_ts, main = "Seasonal Plot", ylab = "", yaxt = "n",
           col = rainbow(8), year.labels = FALSE)
axis_millions()

# ACF after first differencing.
Acf(diff(ampang_ts), lag.max = 24, main = "ACF (first difference)")

# PACF after first differencing.
Pacf(diff(ampang_ts), lag.max = 24, main = "PACF (first difference)")

dev.off()

# STL decomposition.
decomp <- stl(ampang_ts, s.window = "periodic", robust = TRUE)
png(fig("eda", "stl_decomposition.png"), width = 800, height = 600, res = 130)
plot(decomp, main = "STL Decomposition")
dev.off()

# Standalone seasonal plot.
png(fig("eda", "seasonal_plot.png"), width = 800, height = 600, res = 130)
par(mar = c(4, 4, 3, 1))
seasonplot(ampang_ts, main = "Seasonal Plot", ylab = "", yaxt = "n",
           col = rainbow(8), year.labels = TRUE)
axis_millions()
dev.off()

# Lag plots for lags 1 through 12.
png(fig("eda", "lag_plot.png"), width = 800, height = 800, res = 130)
par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
for (k in 1:12) {
  d <- ts.intersect(y = ampang_ts, y_lag = stats::lag(ampang_ts, -k))
  plot(d[, "y_lag"], d[, "y"], xlab = paste0("lag ", k), ylab = "",
       main = paste0("lag ", k), pch = 20, col = "steelblue")
}
dev.off()

# Seasonal strength from the STL decomposition.
seasonal_strength <- max(0, 1 - var(decomp$time.series[, "remainder"]) /
                            var(decomp$time.series[, "seasonal"] + decomp$time.series[, "remainder"]))
cat("Seasonal strength:", round(seasonal_strength, 3), "\n")

# Stationarity tests.
cat("\n--- Level series ---\n")
print(adf.test(ampang_ts))
print(kpss.test(ampang_ts))

cat("\n--- First difference ---\n")
print(adf.test(diff(ampang_ts)))
print(kpss.test(diff(ampang_ts)))