# Replace the MCO-period collapse with a trend and seasonal estimate.
# Keep the full series and use the known disruption window.

source("scripts/setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full.rds")
monthly_df <- readRDS("data/ampang_monthly_full_df.rds")

mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly_df$month >= mco_start & monthly_df$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

# Fill the gap temporarily so STL can estimate seasonality.
temp_filled <- monthly_df$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)
stl_temp <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

# Interpolate the trend between seasonally adjusted endpoint values.
i_pre  <- which(monthly_df$month == mco_start - months(1))
i_post <- which(monthly_df$month == mco_end + months(1))
last_pre_raw   <- monthly_df$rail_lrt_ampang[i_pre]
first_post_raw <- monthly_df$rail_lrt_ampang[i_post]
last_pre   <- last_pre_raw   - seasonal_est[i_pre]
first_post <- first_post_raw - seasonal_est[i_post]
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

cat("Bridging trend from", format(mco_start - months(1)),
    "( raw", last_pre_raw, "-> seasonally adjusted", round(last_pre), ") to",
    format(mco_end + months(1)),
    "( raw", first_post_raw, "-> seasonally adjusted", round(first_post), ")\n")

resolved <- monthly_df$rail_lrt_ampang
resolved[mco_mask] <- bridge_trend + seasonal_est[mco_mask]

comparison <- data.frame(
  month = monthly_df$month[mco_mask],
  original = monthly_df$rail_lrt_ampang[mco_mask],
  resolved = round(resolved[mco_mask])
)
print(comparison)

ampang_ts_resolved <- ts(resolved, start = start(ampang_ts), frequency = 12)
saveRDS(ampang_ts_resolved, "data/ampang_monthly_full_resolved.rds")

cat("\nSaved data/ampang_monthly_full_resolved.rds - length:", length(ampang_ts_resolved),
    "(all", length(ampang_ts), "months retained, MCO window resolved not removed)\n")

# Save a comparison plot.
p_mco <- autoplot(cbind(Original = ampang_ts, Resolved = ampang_ts_resolved)) +
  ggtitle("LRT Ampang: Original vs. MCO-Resolved Series") +
  ylab("Monthly ridership") + scale_y_millions()
print(p_mco)
ggsave(fig("eda", "mco_resolution.png"), p_mco, width = 9, height = 5, dpi = 150)
