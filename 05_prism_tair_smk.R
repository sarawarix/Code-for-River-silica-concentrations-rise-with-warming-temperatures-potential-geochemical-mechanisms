library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
map  <- purrr::map
walk <- purrr::walk


PRISM_DIR <- "/Users/u1322101/data/prism"  
SMK_DIR   <- "data/smk"
WS_RDS    <- "/Users/u1322101/data/smk/ws_all_combined.rds"  

dir.create(PRISM_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SMK_DIR,   recursive = TRUE, showWarnings = FALSE)

YEAR_START <- 1999
YEAR_END   <- 2025  #

MIN_YEARS   <- 3L
MIN_SEASONS <- 4L
MIN_MONTHS  <- 6L

FULL_START <- as.Date("1999-10-01")
FULL_END   <- as.Date("2024-09-30")

water_year <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- as.integer(format(date, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}

wy_period <- function(start_wy, end_wy) {
  list(
    start = as.Date(sprintf("%d-10-01", start_wy - 1L)),
    end   = as.Date(sprintf("%d-09-30", end_wy))
  )
}

all_periods <- c(
  list(list(start = FULL_START, end = FULL_END, label = "WY2000-WY2024")),
  map(2000:2015, function(wy) {
    p <- wy_period(wy, wy + 9L)
    p$label <- sprintf("WY%d-WY%d", wy, wy + 9L)
    p
  })
)


prism_base <- "https://services.nacse.org/prism/data/public/4km"

download_prism_month <- function(year, month, out_dir) {
  ym  <- sprintf("%04d%02d", year, month)
  out <- file.path(out_dir, sprintf("PRISM_tmean_%s.zip", ym))
  if (file.exists(out) && file.size(out) > 10000) return(invisible(out))

  url  <- sprintf("%s/tmean/%s", prism_base, ym)
  resp <- tryCatch(
    download.file(url, destfile = out, mode = "wb", quiet = TRUE),
    error = function(e) -1L
  )

  if (resp == 0L && file.exists(out) && file.size(out) > 10000) {
    cat(sprintf("  Downloaded %s-%02d\n", year, month))
    return(invisible(out))
  }

  warning(sprintf("Failed to download PRISM tmean for %s-%02d", year, month))
  return(invisible(NULL))
}

ym_grid <- expand.grid(year = YEAR_START:YEAR_END, month = 1:12) %>%
  filter(!(year == YEAR_END & month > month(Sys.Date()))) %>%
  filter(!(year == YEAR_START & month < 10)) %>%   
  arrange(year, month)

cat("Months to download:", nrow(ym_grid), "\n")

walk2(ym_grid$year, ym_grid$month, function(yr, mo) {
  download_prism_month(yr, mo, PRISM_DIR)
  Sys.sleep(0.3)  
})


smk_sites <- readRDS(file.path(SMK_DIR, "smk_full.rds")) %>%
  distinct(site_no) %>%
  pull(site_no)
cat("SMK sites to extract:", length(smk_sites), "\n")

ws_all <- if (exists("ws_all_combined")) ws_all_combined else readRDS(WS_RDS)

ws <- ws_all %>%
  mutate(GAGE_ID = str_pad(as.character(GAGE_ID), 8, pad = "0")) %>%
  filter(GAGE_ID %in% smk_sites) %>%
  st_transform(4326)   

zip_files <- list.files(PRISM_DIR, pattern = "\\.zip$", full.names = TRUE)
zip_files <- zip_files[file.size(zip_files) > 10000]
cat("PRISM zip files:   ", length(zip_files), "\n\n")

prism_monthly <- map_dfr(zip_files, function(zf) {
  ym <- str_extract(basename(zf), "\\d{6}")
  yr <- as.integer(substr(ym, 1, 4))
  mo <- as.integer(substr(ym, 5, 6))
  if (mo == 1) cat("  Extracting year:", yr, "\n")

  td  <- tempdir()
  tif <- tryCatch({
    files <- unzip(zf, exdir = td, overwrite = TRUE)
    hits  <- files[grepl("\\.tif$", files, ignore.case = TRUE)]
    if (length(hits) == 0) NA_character_ else hits[1]
  }, error = function(e) NA_character_)

  if (is.na(tif) || !file.exists(tif)) return(NULL)

  r <- tryCatch(terra::rast(tif), error = function(e) NULL)
  if (is.null(r)) return(NULL)

  ws_r <- st_transform(ws, crs = terra::crs(r))
  vals <- exact_extract(r, ws_r, fun = "mean", progress = FALSE)

  tibble(
    site_no     = ws$GAGE_ID,
    year        = yr,
    month       = mo,
    tmean_C     = vals,
    provisional = yr >= 2024,
    Date        = as.Date(sprintf("%04d-%02d-15", yr, mo))
  )
}) %>%
  filter(!is.na(tmean_C))



saveRDS(prism_monthly, file.path(SMK_DIR, "prism_monthly.rds"))


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

run_smk_tair <- function(df_window) {
  empty <- function(reason) tibble(
    tau = NA_real_, p_value = NA_real_, slope = NA_real_,
    method = "rkt_seasonal_median", reason = reason,
    n_years = 0L, n_seasons = 0L, n_months = 0L
  )

  if (!nrow(df_window) || !any(is.finite(df_window$tmean_C)))
    return(empty("no_data"))

  mdf <- df_window %>%
    filter(is.finite(tmean_C)) %>%
    group_by(year, month) %>%
    summarise(val = median(tmean_C, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(val)) %>%
    arrange(year, month)

  n_years   <- n_distinct(mdf$year)
  n_seasons <- n_distinct(mdf$month)
  n_months  <- nrow(mdf)

  if (n_years   < MIN_YEARS)   return(empty(sprintf("too_few_years(%d)",   n_years)))
  if (n_seasons < MIN_SEASONS) return(empty(sprintf("too_few_seasons(%d)", n_seasons)))
  if (n_months  < MIN_MONTHS)  return(empty(sprintf("too_few_months(%d)",  n_months)))
  if (length(unique(mdf$val)) < 2L) return(empty("no_variation"))

  tnum <- as.numeric(mdf$year) + (as.numeric(mdf$month) - 0.5) / 12

  res <- tryCatch(
    .rkt_portable(tnum = tnum, y = mdf$val, block = mdf$month, correct = TRUE),
    error = function(e) e
  )
  if (inherits(res, "error")) return(tibble(
    tau = NA_real_, p_value = NA_real_, slope = NA_real_,
    method   = "rkt_seasonal_median",
    reason   = paste0("rkt_error: ", conditionMessage(res)),
    n_years  = n_years, n_seasons = n_seasons, n_months = n_months
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

sites <- sort(unique(prism_monthly$site_no))
cat("Running SMK on tmean for", length(sites), "sites x",
    length(all_periods), "periods...\n")

smk_tair <- map_dfr(seq_along(sites), function(i) {
  if (i %% 100 == 0) cat("  Site", i, "of", length(sites), "\n")
  sdat <- prism_monthly %>% filter(site_no == sites[i])

  map_dfr(all_periods, function(p) {
    win <- sdat %>% filter(Date >= p$start, Date <= p$end)
    bind_cols(
      tibble(
        site_no      = sites[i],
        variable     = "tmean_C",
        period_start = p$start,
        period_end   = p$end,
        period_label = p$label
      ),
      run_smk_tair(win)
    )
  })
}) %>%
  mutate(
    wy_start  = water_year(period_start + 1),
    wy_end    = water_year(period_end),
    signif_05 = !is.na(p_value) & p_value <= 0.05,
    direction = case_when(
      is.na(tau) ~ NA_character_,
      tau > 0    ~ "increasing",
      tau < 0    ~ "decreasing",
      TRUE       ~ "no trend"
    )
  ) %>%
  arrange(site_no, period_start)

saveRDS(smk_tair, file.path(SMK_DIR, "smk_tair.rds"))

smk_tair %>%
  filter(period_label == "WY2000-WY2024") %>%
  summarise(
    n_sites     = n_distinct(site_no),
    n_valid     = sum(!is.na(tau)),
    pct_warming = round(100 * mean(tau > 0, na.rm = TRUE), 1),
    med_tau     = round(median(tau, na.rm = TRUE), 3),
    med_slope   = round(median(slope, na.rm = TRUE), 4)
  ) %>%
  print()


smk_geo <- readRDS(file.path(SMK_DIR, "smk_geo.rds"))

smk_geo_tair <- bind_rows(
  smk_geo %>% mutate(across(everything(), as.character)),
  smk_tair %>% mutate(across(everything(), as.character))
) %>%
  mutate(
    tau       = as.numeric(tau),
    p_value   = as.numeric(p_value),
    slope     = as.numeric(slope),
    signif_05 = as.logical(signif_05)
  )

saveRDS(smk_geo_tair, file.path(SMK_DIR, "smk_geo_tair.rds"))

