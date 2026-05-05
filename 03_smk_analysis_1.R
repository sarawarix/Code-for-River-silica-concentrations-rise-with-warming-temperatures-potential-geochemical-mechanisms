library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
WQP_FILE <- "data/wqp_download/wqp_screened_all.rds"
CQ_FILE  <- "data/wqp_download/cq_master.rds"
OUT_DIR  <- "data/smk"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MW_SiO2 <- 60.0843
MW_Si   <- 28.0855

MIN_YEARS   <- 3L
MIN_SEASONS <- 4L
MIN_MONTHS  <- 6L

water_year <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- as.integer(format(date, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

wy_period <- function(start_wy, end_wy) {
  c(
    as.character(as.Date(sprintf("%d-10-01", start_wy - 1L))),
    as.character(as.Date(sprintf("%d-09-30", end_wy)))
  )
}


wqp_raw <- readRDS(WQP_FILE) %>%
  mutate(
    Date           = as.Date(Activity_StartDate),
    site_no        = str_remove(Location_Identifier, "^USGS-"),
    Result_Measure = suppressWarnings(as.numeric(Result_Measure)),
    char_lower     = tolower(Result_Characteristic)
  ) %>%
  filter(
    !is.na(Date),
    !is.na(Result_Measure),
    Result_Measure >= 0,
    Date >= as.Date("2000-01-01"),
    Date <= as.Date("2025-09-30")
  )

stopifnot(nrow(wqp_raw) > 0)

wqp_long <- wqp_raw %>%
  mutate(
    variable = case_when(
      char_lower %in% c("silica", "silicon") ~ "Si_mgL",
      char_lower == "ph"                     ~ "pH",
      TRUE                                   ~ NA_character_
    ),
    value = case_when(
      char_lower == "silica"  & Result_MeasureUnit == "mg/L"                 ~ Result_Measure,
      char_lower == "silicon" & Result_MeasureUnit == "ug/L"                 ~
        Result_Measure / 1000 * (MW_SiO2 / MW_Si),
      char_lower %in% c("calcium", "sodium") & Result_MeasureUnit == "mg/L" ~ Result_Measure,
      char_lower == "ph"      & Result_MeasureUnit == "standard units"       ~ Result_Measure,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(variable), !is.na(value), is.finite(value)) %>%
  filter(!(variable == "Si_mgL" & (value <= 0 | value > 50))) %>%
  select(site_no, Date, variable, value)

stopifnot(nrow(wqp_long) > 0)


cq <- readRDS(CQ_FILE) %>%
  mutate(
    Date    = as.Date(Date),
    site_no = str_pad(as.character(site_no), 8, pad = "0")
  ) %>%
  filter(Date >= as.Date("2000-01-01"), Date <= as.Date("2025-09-30"))
stopifnot(nrow(cq) > 0)

cq_long <- cq %>%
  select(site_no, Date, Si_flux_mgs, Flow, Wtemp) %>%
  rename(Q_cfs = Flow, Wtemp_C = Wtemp) %>%
  pivot_longer(
    cols      = c(Si_flux_mgs, Q_cfs, Wtemp_C),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value), is.finite(value), value > 0)


all_long <- bind_rows(wqp_long, cq_long)

stopifnot(nrow(all_long) > 0)

.rkt_portable <- function(tnum, y, block, correct = TRUE) {
  fun <- get("rkt", envir = asNamespace("rkt"))
  fml <- names(formals(fun))
  args <- list(y = y, block = as.integer(block), correct = correct)
  if      ("date" %in% fml) args <- c(list(date = tnum), args)
  else if ("x"    %in% fml) args <- c(list(x    = tnum), args)
  else if ("z"    %in% fml) args <- c(list(z    = tnum), args)
  else stop("rkt::rkt has no recognised time argument (date/x/z)")
  do.call(fun, args)
}

run_smk <- function(df_window) {

  empty_row <- function(reason) {
    tibble(tau = NA_real_, p_value = NA_real_, slope = NA_real_,
           method = "rkt_seasonal_median", reason = reason,
           n_years = 0L, n_seasons = 0L, n_months = 0L)
  }

  if (nrow(df_window) == 0 || !any(is.finite(df_window$value)))
    return(empty_row("no_data"))

  mdf <- df_window %>%
    mutate(year = year(Date), month = month(Date)) %>%
    filter(is.finite(value)) %>%
    group_by(year, month) %>%
    summarise(val = median(value, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(val)) %>%
    arrange(year, month)

  if (nrow(mdf) == 0) return(empty_row("empty_after_monthly_agg"))

  n_years   <- n_distinct(mdf$year)
  n_seasons <- n_distinct(mdf$month)
  n_months  <- nrow(mdf)

  if (n_years   < MIN_YEARS)   return(empty_row(sprintf("too_few_years(%d)",   n_years)))
  if (n_seasons < MIN_SEASONS) return(empty_row(sprintf("too_few_seasons(%d)", n_seasons)))
  if (n_months  < MIN_MONTHS)  return(empty_row(sprintf("too_few_months(%d)",  n_months)))
  if (length(unique(mdf$val)) < 2L) return(empty_row("no_variation"))

  tnum <- as.numeric(mdf$year) + (as.numeric(mdf$month) - 0.5) / 12

  res <- tryCatch(
    .rkt_portable(tnum = tnum, y = mdf$val, block = mdf$month, correct = TRUE),
    error = function(e) e
  )

  if (inherits(res, "error")) return(tibble(
    tau = NA_real_, p_value = NA_real_, slope = NA_real_,
    method  = "rkt_seasonal_median",
    reason  = paste0("rkt_error: ", conditionMessage(res)),
    n_years = n_years, n_seasons = n_seasons, n_months = n_months
  ))

  slope <- unname(
    if      ("B"    %in% names(res)) res$B
    else if ("beta" %in% names(res)) res$beta
    else NA_real_
  )

  tibble(tau = unname(res$tau), p_value = unname(res$sl), slope = slope,
         method = "rkt_seasonal_median", reason = NA_character_,
         n_years = n_years, n_seasons = n_seasons, n_months = n_months)
}

run_smk_all <- function(data_long, periods) {

  sites     <- sort(unique(data_long$site_no))
  var_names <- sort(unique(data_long$variable))

  cat("  Sites:    ", length(sites), "\n")
  cat("  Variables:", paste(var_names, collapse = ", "), "\n")
  cat("  Periods:  ", length(periods), "\n\n")

  results <- map_dfr(seq_along(sites), function(i) {
    if (i %% 100 == 0) cat("  Site", i, "of", length(sites), "\n")

    s         <- sites[i]
    site_data <- data_long %>% filter(site_no == s)

    map_dfr(var_names, function(v) {
      vdata <- site_data %>% filter(variable == v)
      if (nrow(vdata) == 0) return(NULL)

      map_dfr(periods, function(p) {
        start <- as.Date(p[1])
        end   <- as.Date(p[2])
        win   <- vdata %>% filter(Date >= start, Date <= end)
        bind_cols(
          tibble(site_no = s, variable = v,
                 period_start = start, period_end = end),
          run_smk(win)
        )
      })
    })
  })

  if (nrow(results) == 0) {
    warning("run_smk_all: no results returned — check that all_long is non-empty and date window overlaps data")
    return(results)
  }

  results %>%
    mutate(
      signif_05 = !is.na(p_value) & p_value <= 0.05,
      direction = case_when(
        is.na(tau) ~ NA_character_,
        tau > 0    ~ "increasing",
        tau < 0    ~ "decreasing",
        TRUE       ~ "no trend"
      )
    ) %>%
    arrange(variable, site_no, period_start)
}


full_period <- list(wy_period(2000L, 2025L))

smk_full <- run_smk_all(all_long, full_period) %>%
  mutate(period_label = "WY2000-WY2025")

smk_full %>%
  filter(variable == "Si_mgL", !is.na(tau)) %>%
  summarise(
    n_sites          = n(),
    pct_increasing   = round(100 * mean(tau > 0, na.rm = TRUE), 1),
    pct_sig_increase = round(100 * mean(signif_05 & tau > 0, na.rm = TRUE), 1),
    pct_sig_decrease = round(100 * mean(signif_05 & tau < 0, na.rm = TRUE), 1),
    median_tau       = round(median(tau, na.rm = TRUE), 4),
    median_slope     = round(median(slope, na.rm = TRUE), 4),
    q25_slope        = round(quantile(slope, 0.25, na.rm = TRUE), 4),
    q75_slope        = round(quantile(slope, 0.75, na.rm = TRUE), 4)
  ) %>%
  print()

sig_increase_full <- smk_full %>%
  filter(variable == "Si_mgL", !is.na(tau), signif_05, tau > 0)


saveRDS(smk_full, file.path(OUT_DIR, "smk_full.rds"))
write.csv(smk_full, file.path(OUT_DIR, "smk_full.csv"), row.names = FALSE)

# ROLLING 10-WATER-YEAR WINDOWS
window_start_wys <- 2000:2015

periods_rolling <- purrr::map(window_start_wys, function(wy_start) {
  wy_period(wy_start, wy_start + 9L)
})

purrr::walk2(window_start_wys, periods_rolling, function(wy, p) {
  cat(sprintf("  WY%d-WY%d  (%s to %s)\n", wy, wy + 9L, p[1], p[2]))
})
cat("\n")

smk_rolling <- run_smk_all(all_long, periods_rolling) %>%
  mutate(
    wy_start     = water_year(period_start + 1),
    wy_end       = water_year(period_end),
    period_label = sprintf("WY%d-WY%d", wy_start, wy_end),
    mid_wy       = wy_start + 5L
  )

si_summary <- smk_rolling %>%
  filter(variable == "Si_mgL", !is.na(tau)) %>%
  group_by(period_label, wy_start, wy_end, mid_wy) %>%
  summarise(
    n_sites          = n(),
    median_tau       = round(median(tau, na.rm = TRUE), 4),
    mean_tau         = round(mean(tau,   na.rm = TRUE), 4),
    pct_increasing   = round(100 * mean(tau > 0,                    na.rm = TRUE), 1),
    pct_sig_increase = round(100 * mean(signif_05 & tau > 0,        na.rm = TRUE), 1),
    pct_sig_decrease = round(100 * mean(signif_05 & tau < 0,        na.rm = TRUE), 1),
    consistency_idx  = round(pct_sig_increase - pct_sig_decrease, 1),
    .groups          = "drop"
  ) %>%
  arrange(desc(consistency_idx))

print(si_summary, n = 20)
peak_window <- si_summary %>% slice_head(n = 1)


sig_peak <- smk_rolling %>%
  filter(variable == "Si_mgL",
         period_label == peak_window$period_label,
         !is.na(tau), signif_05, tau > 0)

si_tau_plot <- smk_rolling %>%
  filter(variable == "Si_mgL", !is.na(tau)) %>%
  mutate(period_label = factor(period_label,
                               levels = unique(period_label[order(wy_start)])))

p_violin <- ggplot(si_tau_plot, aes(x = period_label, y = tau, fill = mid_wy)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(alpha = 0.7, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.4,
               colour = "grey20", fill = "white", alpha = 0.6) +
  scale_fill_gradient(low = "#4393c3", high = "#d6604d",
                      name = "Window\nmidpoint\n(water year)") +
  labs(x = "10-water-year window", y = "Kendall \u03c4",
       title = "Si_mgL SMK \u03c4, rolling 10-water-year windows") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(OUT_DIR, "smk_tau_distributions.png"),
       p_violin, width = 14, height = 6, dpi = 150)

p_consistency <- ggplot(si_summary, aes(x = wy_start, y = consistency_idx)) +
  geom_line(colour = "#2166ac", linewidth = 1) +
  geom_point(aes(colour = consistency_idx), size = 3) +
  scale_colour_gradient2(low = "#d6604d", mid = "grey80", high = "#2166ac",
                         midpoint = 0, name = "Consistency\nindex") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(breaks = window_start_wys, name = "Window start (water year)") +
  labs(y = "Consistency index (pp)",
       title = "Si_mgL trend consistency, rolling 10-water-year windows") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "smk_consistency_index.png"),
       p_consistency, width = 10, height = 5, dpi = 150)

saveRDS(smk_rolling, file.path(OUT_DIR, "smk_rolling.rds"))
write.csv(smk_rolling, file.path(OUT_DIR, "smk_rolling.csv"),        row.names = FALSE)
write.csv(si_summary,  file.path(OUT_DIR, "smk_decade_summary.csv"), row.names = FALSE)

wqp_si_full <- all_long %>%
  filter(variable == "Si_mgL",
         Date >= as.Date("1999-10-01"),
         Date <= as.Date("2025-09-30"))


run_smk_season <- function(df) {
  empty <- tibble(tau = NA_real_, p_value = NA_real_, slope = NA_real_,
                  n_years = 0L, n_months = 0L)
  if (nrow(df) == 0) return(empty)
  mdf <- df %>%
    mutate(year = year(Date), month = month(Date)) %>%
    filter(is.finite(value)) %>%
    group_by(year, month) %>%
    summarise(val = median(value, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(val)) %>%
    arrange(year, month)
  n_years  <- n_distinct(mdf$year)
  n_months <- nrow(mdf)
  if (n_years < MIN_YEARS || n_months < MIN_MONTHS) return(empty)
  if (length(unique(mdf$val)) < 2) return(empty)
  tnum <- as.numeric(mdf$year) + (as.numeric(mdf$month) - 0.5) / 12
  res <- tryCatch(
    .rkt_portable(tnum = tnum, y = mdf$val, block = mdf$month, correct = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) return(empty)
  slope <- unname(
    if      ("B"    %in% names(res)) res$B
    else if ("beta" %in% names(res)) res$beta
    else NA_real_
  )
  tibble(tau = unname(res$tau), p_value = unname(res$sl), slope = slope,
         n_years = n_years, n_months = n_months)
}

sites_seasonal <- sort(unique(wqp_si_full$site_no))

smk_seasonal <- map_dfr(sites_seasonal, function(s) {
  df <- wqp_si_full %>% filter(site_no == s)
  bind_rows(
    run_smk_season(df %>% filter(month(Date) %in% 6:8)) %>%
      mutate(site_no = s, season = "Summer (JJA)"),
    run_smk_season(df %>% filter(month(Date) %in% c(12, 1, 2))) %>%
      mutate(site_no = s, season = "Winter (DJF)")
  )
}) %>%
  filter(!is.na(tau)) %>%
  mutate(signif_05 = p_value <= 0.05)


seasonal_wide <- smk_seasonal %>%
  select(site_no, season, tau) %>%
  pivot_wider(names_from = season, values_from = tau) %>%
  filter(!is.na(`Summer (JJA)`), !is.na(`Winter (DJF)`))

wilcox_overall <- wilcox.test(
  seasonal_wide$`Summer (JJA)`,
  seasonal_wide$`Winter (DJF)`,
  paired = TRUE, alternative = "greater"
)

saveRDS(smk_seasonal, file.path(OUT_DIR, "smk_seasonal.rds"))
saveRDS(all_long, file.path(OUT_DIR, "all_long.rds"))