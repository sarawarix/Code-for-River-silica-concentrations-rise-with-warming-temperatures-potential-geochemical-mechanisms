library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)
library(stringr)

WQP_FILE   <- "data/wqp_download/wqp_screened_all.rds"
PRISM_FILE <- "data/smk/prism_monthly.rds"
OUT_DIR    <- "data/smk"

MW_SiO2 <- 60.0843
MW_Si   <- 28.0855
MAX_LAG <- 10L
MIN_YRS <- 11L   

ANAL_START <- as.Date("1990-10-01")
ANAL_END   <- as.Date("2025-09-30")

water_year <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- as.integer(format(date, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

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

prism <- readRDS(PRISM_FILE) %>%
  filter(Date >= as.Date("1989-10-01")) %>%   
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

si_lm   <- lm(med_si   ~ wy, data = ts_dat)
tair_lm <- lm(med_tair ~ wy, data = ts_dat)

sites <- intersect(unique(si_ann$site_no), unique(tair_ann$site_no))

all_lag_cors <- map_dfr(sites, function(s) {
  si   <- si_ann   %>% filter(site_no == s)
  tair <- tair_ann %>% filter(site_no == s)

  map_dfr(0:MAX_LAG, function(lag) {
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

saveRDS(all_lag_cors, file.path(OUT_DIR, "lag_analysis.rds"))

optimal_lags <- all_lag_cors %>%
  group_by(site_no) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    direction = if_else(r > 0, "positive", "negative"),
    sig       = p_value <= 0.05
  )

optimal_lags %>%
  filter(direction == "positive") %>%
  count(lag) %>%
  print()

saveRDS(optimal_lags, file.path(OUT_DIR, "optimal_lags.rds"))
