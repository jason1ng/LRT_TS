# Aggregate daily Ampang LRT ridership into monthly totals.
# Source: data.gov.my Daily Public Transport Ridership dataset
# (https://data.gov.my/data-catalogue/ridership_headline)

source("scripts/setup.R")

raw <- read.csv("data/ridership_headline.csv", stringsAsFactors = FALSE)
raw$date <- as.Date(raw$date)

# Sum trips by calendar month.
monthly <- raw %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE),
    n_days = n()
  ) %>%
  ungroup()

# Exclude incomplete months.
monthly <- monthly %>% filter(n_days >= 28)

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

# Create a monthly time series.
start_year  <- year(min(monthly$month))
start_month <- month(min(monthly$month))
ampang_ts <- ts(monthly$rail_lrt_ampang, start = c(start_year, start_month), frequency = 12)

ts.plot(
  ampang_ts / 1000,
  main = "Monthly LRT Ampang Ridership",
  xlab = "Year",
  ylab = "Ridership ('000)",
  lwd = 1
)

# Basic checks.
nrow(monthly)
range(monthly$month)
monthly %>% arrange(month) %>% select(month, n_days) %>% print(n = Inf)

# 2020 minimum.
covid1 <- monthly %>% filter(year(month) == 2020)
covid1[which.min(covid1$rail_lrt_ampang), c("month", "rail_lrt_ampang")]

# 2021 minimum.
covid2 <- monthly %>% filter(year(month) == 2021)
covid2[which.min(covid2$rail_lrt_ampang), c("month", "rail_lrt_ampang")]

# 2022 range.
range2022 <- monthly %>% filter(year(month) == 2022)
summary(range2022$rail_lrt_ampang)
saveRDS(ampang_ts, "data/ampang_monthly_full.rds")
saveRDS(monthly, "data/ampang_monthly_full_df.rds")

cat("Saved data/ampang_monthly_full.rds - length:", length(ampang_ts), "\n")
