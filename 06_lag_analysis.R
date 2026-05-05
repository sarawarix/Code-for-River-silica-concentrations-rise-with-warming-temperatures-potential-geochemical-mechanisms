# =============================================================
# SCRIPT 06 — Si–Temperature Lag Analysis
#
# PURPOSE:
#   For each site with >= 11 years of paired annual silica and
#   air temperature observations, find the lag (0–10 years)
#   between PRISM mean annual temperature and annual median Si
#   that maximises the Pearson correlation |r|.
#
# EXCLUSION:
#   Sites where the optimal correlation is negative are excluded
#   from lag plots (per methods), but retained in the full results.
#
# INPUT:
#   - wqp_screened_all.rds     — WQP data
#   - data/smk/prism_monthly.rds — monthly PRISM tmean per site
#
# OUTPUT:
#   - data/smk/lag_analysis.rds   — one row per site, all lags tested
#   - data/smk/optimal_lags.rds   — one row per site, optimal lag only
#   - data/smk/ts_dat.rds         — annual cross-site medians for Figure 2A
# =============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)
library(stringr)

# -------------------------------------------------------------
# 0. CONFIG
# -------------------------------------------------------------
WQP_FILE   <- "data/wqp_download/wqp_screened_all.rds"
PRISM_FILE <- "data/smk/prism_monthly.rds"
OUT_DIR    <- "data/smk"

MW_SiO2 <- 60.0843
MW_Si   <- 28.0855
MAX_LAG <- 10L
MIN_YRS <- 11L   # per methods: "at least 11 continuous years"

# Analysis window: WY2000–WY2024
ANAL_START <- as.Date("1990-10-01")
ANAL_END   <- as.Date("2025-09-30")

water_year <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- as.integer(format(date, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

# -------------------------------------------------------------
# 1. ANNUAL MEDIAN SILICA PER SITE
# Water-year medians; require >= 2 observations per water year
# -------------------------------------------------------------
cat("Computing annual median silica...\n")

si_ann <- readRDS(WQP_FILE) %>%
  mutate(
    Date           = as.Date(Activity_StartDate),   # character, direct parse
    site_no        = str_remove(Location_Identifier, "^USGS-"),
    Result_Measure = suppressWarnings(as.numeric(Result_Measure)),
    char_lower     = tolower(Result_Characteristic)
  ) %>%
  filter(
    !is.na(Date), !is.na(Result_Measure),
    char_lower %in% c("silica", "silicon"),
    Date >= ANAL_START, Date <= ANAL_END
  ) %>%
  mutate(
    si_mgL = case_when(
      char_lower == "silica"  & Result_MeasureUnit == "mg/L" ~ Result_Measure,
      char_lower == "silicon" & Result_MeasureUnit == "ug/L" ~
        Result_Measure / 1000 * (MW_SiO2 / MW_Si),
      TRUE ~ NA_real_
    ),
    wy = water_year(Date)
  ) %>%
  filter(!is.na(si_mgL), si_mgL > 0, si_mgL <= 50) %>%
  group_by(site_no, wy) %>%
  summarise(
    si_med = median(si_mgL, na.rm = TRUE),
    n_obs  = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 2)

cat("si_ann rows:", nrow(si_ann), "\n")
cat("si_ann sites:", n_distinct(si_ann$site_no), "\n")
cat("si_ann wy range:", range(si_ann$wy), "\n")

cat("Annual Si site-years:", nrow(si_ann), "\n")
cat("Unique sites:        ", n_distinct(si_ann$site_no), "\n\n")

# -------------------------------------------------------------
# 2. ANNUAL MEAN TEMPERATURE PER SITE (water-year)
# Also compute site baseline (mean over full period) for anomalies
# -------------------------------------------------------------
cat("Computing annual mean air temperature...\n")

prism <- readRDS(PRISM_FILE) %>%
  filter(Date >= as.Date("1989-10-01")) %>%   # <-- add %>% here
  mutate(wy = water_year(Date))

tair_ann <- prism %>%
  group_by(site_no, wy) %>%
  summarise(tmean = mean(tmean_C, na.rm = TRUE), .groups = "drop")

tair_base <- tair_ann %>%
  group_by(site_no) %>%
  summarise(baseline = mean(tmean, na.rm = TRUE), .groups = "drop")

tair_ann <- tair_ann %>%
  left_join(tair_base, by = "site_no") %>%
  mutate(tair_anom = tmean - baseline)

cat("Annual temp site-years:", nrow(tair_ann), "\n\n")

# -------------------------------------------------------------
# 3. CROSS-SITE ANNUAL MEDIANS FOR FIGURE 2A
# Annual median across all sites with data in that water year
# -------------------------------------------------------------
ts_dat <- si_ann %>%
  inner_join(tair_ann %>% select(site_no, wy, tmean), by = c("site_no", "wy")) %>%
  group_by(wy) %>%
  summarise(
    med_si   = median(si_med,  na.rm = TRUE),
    med_tair = median(tmean,   na.rm = TRUE),
    n_sites  = n(),
    .groups  = "drop"
  )

saveRDS(ts_dat, file.path(OUT_DIR, "ts_dat.rds"))
cat("Saved ts_dat.rds (", nrow(ts_dat), "water years)\n\n")

# Long-term regression stats for Figure 2A annotation
si_lm   <- lm(med_si   ~ wy, data = ts_dat)
tair_lm <- lm(med_tair ~ wy, data = ts_dat)
cat("Si trend:   slope =", round(coef(si_lm)["wy"], 4),
    "mg/L/yr,  p =", round(summary(si_lm)$coefficients["wy", "Pr(>|t|)"], 4), "\n")
cat("Tair trend: slope =", round(coef(tair_lm)["wy"], 4),
    "deg C/yr, p =", round(summary(tair_lm)$coefficients["wy", "Pr(>|t|)"], 4), "\n\n")

# -------------------------------------------------------------
# 4. OPTIMAL LAG ANALYSIS
# For each site x lag (0–10 yrs): correlate annual Si median
# against tair anomaly lagged by that many years.
# Sites need MIN_YRS paired observations at a given lag.
# -------------------------------------------------------------
cat("Running optimal lag analysis (lags 0-", MAX_LAG, ")...\n")

sites <- intersect(unique(si_ann$site_no), unique(tair_ann$site_no))
cat("Sites with both Si and tair:", length(sites), "\n\n")

all_lag_cors <- map_dfr(sites, function(s) {
  si   <- si_ann   %>% filter(site_no == s)
  tair <- tair_ann %>% filter(site_no == s)

  map_dfr(0:MAX_LAG, function(lag) {
    # Si in water year WY paired with tair from WY - lag
    joined <- si %>%
      mutate(tair_wy = wy - lag) %>%
      left_join(
        tair %>% select(site_no, wy, tair_anom),
        by = c("site_no", "tair_wy" = "wy")
      ) %>%
      filter(!is.na(tair_anom), !is.na(si_med), is.finite(si_med))

    if (nrow(joined) < MIN_YRS) return(NULL)

    ct <- tryCatch(
      cor.test(joined$si_med, joined$tair_anom, method = "pearson"),
      error = function(e) NULL
    )
    if (is.null(ct)) return(NULL)

    tibble(
      site_no = s,
      lag     = lag,
      r       = unname(ct$estimate),
      p_value = ct$p.value,
      n_yrs   = nrow(joined)
    )
  })
})

cat("Lag-correlation rows computed:", nrow(all_lag_cors), "\n")
saveRDS(all_lag_cors, file.path(OUT_DIR, "lag_analysis.rds"))

# -------------------------------------------------------------
# 5. OPTIMAL LAG PER SITE
# Select the lag with the highest |r| for each site.
# If multiple lags tie, take the smallest lag.
# -------------------------------------------------------------
optimal_lags <- all_lag_cors %>%
  group_by(site_no) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    direction = if_else(r > 0, "positive", "negative"),
    sig       = p_value <= 0.05
  )

cat("\nOptimal lag summary:\n")
cat("  Sites with optimal lag:    ", nrow(optimal_lags), "\n")
cat("  Positive correlation:      ", sum(optimal_lags$direction == "positive"), "\n")
cat("  Negative correlation:      ", sum(optimal_lags$direction == "negative"), "\n")
cat("  Significant (p<=0.05):     ", sum(optimal_lags$sig), "\n")

cat("\nLag distribution (positive sites):\n")
optimal_lags %>%
  filter(direction == "positive") %>%
  count(lag) %>%
  print()

saveRDS(optimal_lags, file.path(OUT_DIR, "optimal_lags.rds"))
cat("Saved: lag_analysis.rds, optimal_lags.rds\n")
