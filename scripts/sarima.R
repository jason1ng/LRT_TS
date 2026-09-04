# Fit and evaluate a SARIMA model.

source("scripts/setup.R")

y  <- load_series()
sp <- split_series(y)
train <- sp$train
test  <- sp$test

# Stage 1: identification.

# Select seasonal and non-seasonal differencing orders.
D <- nsdiffs(train)
d <- ndiffs(if (D > 0) diff(train, lag = 12) else train, test = "kpss")
cat("Identification: d =", d, ", D =", D, "(seasonal period = 12)\n")

# Inspect the stationary series before fitting candidate models.
stationary_train <- diff(train, differences = d)
if (D > 0) stationary_train <- diff(stationary_train, lag = 12)

png(fig("models", "sarima_identification_acf_pacf.png"), width = 800, height = 400, res = 130)
par(mfrow = c(1, 2))
Acf(stationary_train,  lag.max = 24, main = "ACF (stationary series)")
Pacf(stationary_train, lag.max = 24, main = "PACF (stationary series)")
dev.off()

# Stage 2: estimation and model selection.

# Compare a small grid with fixed differencing orders.
grid <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
results <- data.frame()

for (i in seq_len(nrow(grid))) {
  fit <- tryCatch(
    Arima(train,
          order    = c(grid$p[i], d, grid$q[i]),
          seasonal = list(order = c(grid$P[i], D, grid$Q[i]), period = 12)),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    results <- rbind(results, data.frame(
      p = grid$p[i], d = d, q = grid$q[i],
      P = grid$P[i], D = D, Q = grid$Q[i],
      AICc = fit$aicc, BIC = fit$bic))
  }
}

results <- results[order(results$AICc), ]
cat("\nTop 8 candidate models by AICc:\n")
print(head(results, 8), row.names = FALSE)
write.csv(results, tbl("sarima_grid.csv"), row.names = FALSE)

# Check the AICc gap to the next candidate.
best  <- results[1, ]
label <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[12]",
                  best$p, best$d, best$q, best$P, best$D, best$Q)
cat("\nSelected:", label, "| AICc =", round(best$AICc, 2), "\n")

if (nrow(results) > 1) {
  gap <- results$AICc[2] - results$AICc[1]
  cat("Gap to 2nd-best candidate: AICc +", round(gap, 2),
      if (gap <= 2) "(not a substantial difference - simpler model is defensible)\n"
      else "(substantial - top model clearly preferred)\n")
}

# Refit the selected model.
fit_sarima <- Arima(train,
                     order    = c(best$p, best$d, best$q),
                     seasonal = list(order = c(best$P, best$D, best$Q), period = 12))
print(summary(fit_sarima))
# Approximate coefficient significance using two standard errors.
cat("\nCoefficient significance (|coef| > 2*s.e. approx. p<0.05):\n")
coef_table <- data.frame(
  term = names(coef(fit_sarima)),
  estimate = round(coef(fit_sarima), 4),
  se = round(sqrt(diag(vcov(fit_sarima))), 4)
)
coef_table$significant <- abs(coef_table$estimate) > 2 * coef_table$se
print(coef_table, row.names = FALSE)

# Stage 3: diagnostic checking.

fc_sarima <- forecast(fit_sarima, h = H)

png(fig("models", "sarima_forecast.png"), width = 800, height = 500, res = 130)
plot(fc_sarima, main = paste(label, "- Forecast vs Actual"),
     xlab = "Year", ylab = "Ridership")
lines(test, col = "red", lwd = 2)
legend("topleft", legend = c("Forecast", "Actual (test)"),
       col = c("blue", "red"), lty = 1, lwd = 2, bty = "n")
dev.off()

cat("\nHoldout accuracy:\n")
print(accuracy(fc_sarima, test))

# Test residual autocorrelation with the fitted model degrees of freedom.
cat("\nLjung-Box test (residuals should be white noise, want p > 0.05):\n")
print(lb_test(fit_sarima, 12))
print(lb_test(fit_sarima, LAG_MAX))
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(residuals(fit_sarima), lag.max = 12),
    "/ 12 | lag", LAG_MAX, "=",
    acf_out_of_bounds(residuals(fit_sarima), lag.max = LAG_MAX), "/", LAG_MAX, "\n")

png(fig("models", "sarima_resid.png"), width = 900, height = 700, res = 130)
checkresiduals(fit_sarima)
dev.off()

png(fig("models", "sarima_residuals.png"), width = 900, height = 700, res = 130)
par(mfrow = c(2, 2))
plot(residuals(fit_sarima), main = "Residuals over time", ylab = "")
abline(h = 0, col = "red", lty = 2)
Acf(residuals(fit_sarima), lag.max = LAG_MAX, main = "Residual ACF (should show no spikes)")
hist(residuals(fit_sarima), main = "Residual histogram", xlab = "")
qqnorm(residuals(fit_sarima)); qqline(residuals(fit_sarima), col = "red")
dev.off()

# Compare training and test accuracy.
g <- gap_check(accuracy(fc_sarima, test))
cat("\nOverfitting check (", g$direction, ")\n")
cat("  MASE gap  :", round(g$mase_gap * 100, 1), "% | <= 10%  :", g$within_10pct, "\n")
cat("  RMSE ratio:", round(g$rmse_ratio, 3), "  | <= 1.3x :", g$within_1_3x, "\n")

# Save the model summary.
summ <- model_summary(label, fit_sarima, fc_sarima, test)
print(summ)
write.csv(summ, tbl("summary_sarima.csv"), row.names = FALSE)