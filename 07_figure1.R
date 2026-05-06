library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(patchwork)
library(maps)

BASIN_FILE <- "/Users/u1322101/Library/CloudStorage/OneDrive-UniversityofUtah/Documents_OneDrive/Research Projects/silica_temp/gagesII/Dataset1_BasinID/BasinID.txt"
WQP_FILE   <- "data/wqp_download/wqp_screened_all.rds"
SMK_ROLL   <- "data/smk/smk_rolling.rds"
FIG_DIR    <- "data/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MW_SiO2 <- 60.0843
MW_Si   <- 28.0855

water_year <- function(date) {
  yr <- as.integer(format(date, "%Y"))
  mo <- as.integer(format(date, "%m"))
  ifelse(mo >= 10L, yr + 1L, yr)
}


basins <- read.csv(BASIN_FILE, stringsAsFactors = FALSE, strip.white = TRUE) %>%
  mutate(
    STAID   = str_pad(as.character(STAID), 8, pad = "0"),
    site_no = STAID,
    lat     = LAT_GAGE,
    lon     = LNG_GAGE
  ) %>%
  select(site_no, lat, lon)

smk_rolling <- readRDS(SMK_ROLL)

si_summary <- smk_rolling %>%
  filter(variable == "Si_mgL", !is.na(tau)) %>%
  group_by(period_label, period_start, period_end) %>%
  summarise(
    consistency_idx = mean(
      (!is.na(p_value) & p_value <= 0.05 & tau > 0), na.rm = TRUE
    ) - mean(
      (!is.na(p_value) & p_value <= 0.05 & tau < 0), na.rm = TRUE
    ),
    .groups = "drop"
  )

peak_window <- si_summary %>% slice_max(consistency_idx, n = 1)
PEAK_START  <- peak_window$period_start
PEAK_END    <- peak_window$period_end
PEAK_WY_S   <- water_year(PEAK_START + 1)
PEAK_WY_E   <- water_year(PEAK_END)
PEAK_LABEL  <- sprintf("WY%d\u2013WY%d", PEAK_WY_S, PEAK_WY_E)


wqp_raw <- readRDS(WQP_FILE) %>%
  mutate(
    Date           = as.Date(Activity_StartDate),   
    site_no        = str_remove(Location_Identifier, "^USGS-"),
    Result_Measure = suppressWarnings(as.numeric(Result_Measure)),
    char_lower     = tolower(Result_Characteristic)
  ) %>%
  filter(
    char_lower %in% c("silica", "silicon"),
    !is.na(Date), !is.na(Result_Measure),
    Date >= as.Date("1999-10-01"),   # start of WY2000
    Date <= as.Date("2025-09-30")
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
  filter(!is.na(si_mgL), si_mgL > 0, si_mgL <= 25)


ecdf_dat <- wqp_raw %>%
  group_by(site_no, wy) %>%
  summarise(si_med = median(si_mgL, na.rm = TRUE), .groups = "drop") %>%
  group_by(wy) %>%
  arrange(si_med) %>%
  mutate(ecdf_y = seq_along(si_med) / n()) %>%
  ungroup()

cat("ECDF data: ", n_distinct(ecdf_dat$wy), "water years,",
    n_distinct(ecdf_dat$site_no), "sites\n")


map_dat <- smk_rolling %>%
  filter(
    variable     == "Si_mgL",
    period_start == PEAK_START,
    period_end   == PEAK_END,
    !is.na(tau)
  ) %>%
  left_join(basins, by = "site_no") %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  arrange(abs(tau))   # plot small |tau| first so large values are on top

us_states <- map_data("state")


pA <- ggplot(ecdf_dat,
             aes(x = si_med, y = ecdf_y, group = wy, color = wy)) +

  geom_line(linewidth = 0.4, alpha = 0.80) +

  scale_color_viridis_c(
    option = "plasma",
    name   = "Water year",
    breaks = c(2000, 2005, 2010, 2015, 2020, 2025),
    guide  = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth  = unit(0.35, "cm"),
      barheight = unit(3.5, "cm"),
      direction = "vertical"
    )
  ) +

  scale_x_continuous(limits = c(0, 25),
                     name   = "Site-median annual Si (mg SiO\u2082/L)") +
  scale_y_continuous(name   = "ECDF",
                     breaks = seq(0, 1, 0.25)) +

  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey93"),
    legend.position   = c(0.92, 0.38),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7.5),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA)
  )


pB <- ggplot() +

  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey70", color = "grey40", linewidth = 0.25) +

  geom_point(data = map_dat,
             aes(x = lon, y = lat, color = tau),
             size = 2.5, alpha = 0.88, shape = 16) +

  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "grey92",
    high     = "#B2182B",
    midpoint = 0,
    name     = paste0("SMK \u03c4\n(", PEAK_LABEL, ")"),
    limits   = c(-0.6, 0.6),
    oob      = scales::squish,
    breaks   = c(-0.3, 0.0, 0.3),
    guide    = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth  = unit(0.35, "cm"),
      barheight = unit(3.0, "cm"),
      direction = "vertical"
    )
  ) +

  coord_fixed(1.3, xlim = c(-125, -65), ylim = c(24, 50)) +
  labs(x = NULL, y = NULL) +

  theme_bw(base_size = 10) +
  theme(
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = c(0.92, 0.32),
    legend.direction  = "vertical",
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7.5),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA)
  )


fig1 <- (pA / pB +
           plot_layout(heights = c(1, 1.2))) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag          = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.02, 0.97)
    )
  )

# Save
ggsave(file.path(FIG_DIR, "figure1.pdf"), fig1, width = 7, height = 8)
ggsave(file.path(FIG_DIR, "figure1.png"), fig1, width = 7, height = 8, dpi = 300)

fig1
