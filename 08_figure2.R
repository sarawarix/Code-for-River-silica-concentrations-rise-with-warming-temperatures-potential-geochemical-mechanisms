# =============================================================
# SCRIPT 08 — Figure 2
#
# PANEL A: Median annual Si (blue) and median air temperature
#   (red) across all 191 sites, WY2000–WY2024, with dual y-axes
#   and regression slopes annotated.
#
# PANEL B: Histogram of optimal Si–temperature lag times for
#   sites with positive correlations (Figure 2B in paper).
#
# PANEL C: Violin + boxplot of SMK tau by dominant watershed
#   lithology, full period WY2000–WY2024, with summer/winter
#   breakdown (Figure 2C in paper).
#
# INPUT:
#   - data/smk/ts_dat.rds           — from 06_lag_analysis.R
#   - data/smk/optimal_lags.rds     — from 06_lag_analysis.R
#   - data/smk/smk_geo.rds          — from 04_geology_assignment.R
#   - data/smk/smk_seasonal_geo.rds — from 04_geology_assignment.R
#
# OUTPUT:
#   - data/figures/figure2.pdf  (7 x 10 inches)
#   - data/figures/figure2.png
# =============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)

FIG_DIR <- "data/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Geology color palette (matches paper description)
geo_colors <- c(
  "Igneous"       = "#4A85B8",
  "Metamorphic"   = "#5A9E6F",
  "Sandstone"     = "#8B72B0",
  "Siltstone"     = "#C9B84C",
  "Clay/Mudstone" = "#C07B7B"
)

season_colors <- c(
  "Summer (JJA)" = "#E07B39",
  "Winter (DJF)" = "#4A85B8"
)

# =============================================================
# PANEL A — Si + Temperature time series
# =============================================================
ts_dat <- readRDS("data/smk/ts_dat.rds")

# Restrict to WY2000–WY2024 for plot (match methods window)
ts_plot <- ts_dat %>% filter(wy >= 2000, wy <= 2024)

# Regressions for slope annotations
si_lm   <- lm(med_si   ~ wy, data = ts_plot)
tair_lm <- lm(med_tair ~ wy, data = ts_plot)

si_slope <- coef(si_lm)["wy"]
si_p     <- summary(si_lm)$coefficients["wy", "Pr(>|t|)"]
t_slope  <- coef(tair_lm)["wy"]
t_p      <- summary(tair_lm)$coefficients["wy", "Pr(>|t|)"]

cat("Figure 2A regressions:\n")
cat(sprintf("  Si:   slope = %+.3f mg/L/yr,  p = %.4f\n", si_slope, si_p))
cat(sprintf("  Tair: slope = %+.4f deg/yr, p = %.4f\n",   t_slope,  t_p))

# Dual-axis scaling: map tair onto Si axis
si_min <- min(ts_plot$med_si,   na.rm = TRUE)
si_max <- max(ts_plot$med_si,   na.rm = TRUE)
t_min  <- min(ts_plot$med_tair, na.rm = TRUE)
t_max  <- max(ts_plot$med_tair, na.rm = TRUE)

scale_fac <- (si_max - si_min) / (t_max - t_min)
t_offset  <- mean(ts_plot$med_si, na.rm = TRUE) -
  mean(ts_plot$med_tair, na.rm = TRUE) * scale_fac

ts_plot <- ts_plot %>%
  mutate(tair_scaled = med_tair * scale_fac + t_offset)

y_lo <- min(ts_plot$med_si, ts_plot$tair_scaled, na.rm = TRUE) * 0.97
y_hi <- max(ts_plot$med_si, ts_plot$tair_scaled, na.rm = TRUE) * 1.03

si_label   <- sprintf("Si: %+.3f mg/L yr\u207b\u00b9,  p = %.4f",  si_slope, si_p)
tair_label <- sprintf("T\u2090\u1d35\u1d63: %+.4f \u00b0C yr\u207b\u00b9, p = %.4f", t_slope, t_p)

pA <- ggplot(ts_plot, aes(x = wy)) +
  geom_line(aes(y = tair_scaled, color = "Air temperature"),
            linewidth = 0.55, alpha = 0.45) +
  geom_line(aes(y = med_si, color = "Silica"),
            linewidth = 0.55, alpha = 0.45) +
  geom_point(aes(y = tair_scaled, color = "Air temperature"),
             size = 2.0, alpha = 0.90, shape = 16) +
  geom_point(aes(y = med_si, color = "Silica"),
             size = 2.0, alpha = 0.85, shape = 16) +
  # Regression lines
  geom_smooth(aes(y = tair_scaled, color = "Air temperature"),
              method = "lm", se = FALSE, linewidth = 0.7, linetype = "dashed") +
  geom_smooth(aes(y = med_si, color = "Silica"),
              method = "lm", se = FALSE, linewidth = 0.7, linetype = "dashed") +
  # Annotation
  annotate("text", x = 2024, y = Inf, label = si_label,
           hjust = 1, vjust = 4,   size = 2.8, color = "#2C3E6B") +
  annotate("text", x = 2024, y = Inf, label = tair_label,
           hjust = 1, vjust = 5.8, size = 2.8, color = "#C0392B") +
  scale_color_manual(
    values = c("Silica" = "#2C3E6B", "Air temperature" = "#C0392B"),
    name   = NULL
  ) +
  scale_y_continuous(
    name     = "Median Si (mg SiO\u2082 L\u207b\u00b9)",
    limits   = c(y_lo, y_hi),
    expand   = expansion(mult = 0.05),
    sec.axis = sec_axis(~ (. - t_offset) / scale_fac,
                        name = "Median air temperature (\u00b0C)")
  ) +
  scale_x_continuous(breaks = seq(2000, 2024, by = 4),
                     name   = "Water year") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "grey93"),
    legend.position    = c(0.15, 0.88),
    legend.background  = element_rect(fill = alpha("white", 0.8), colour = NA),
    legend.text        = element_text(size = 8.5),
    axis.title.y.left  = element_text(color = "#2C3E6B", face = "bold", size = 9),
    axis.text.y.left   = element_text(color = "#2C3E6B"),
    axis.title.y.right = element_text(color = "#C0392B", face = "bold", size = 9),
    axis.text.y.right  = element_text(color = "#C0392B"),
    plot.margin        = ggplot2::margin(4, 4, 0, 4)
  )

# =============================================================
# PANEL B — Optimal lag histogram (positive correlation sites)
# =============================================================
optimal_lags <- readRDS("data/smk/optimal_lags.rds")

pos_dat <- optimal_lags %>%
  filter(direction == "positive") %>%
  mutate(
    sig_label = if_else(sig, "p \u2264 0.05", "p > 0.05"),
    sig_label = factor(sig_label, levels = c("p > 0.05", "p \u2264 0.05"))
  )

cat("\nPositive correlation sites:", nrow(pos_dat), "\n")
cat("Significant (p<=0.05):", sum(pos_dat$sig), "\n")

pB <- ggplot(pos_dat, aes(x = factor(lag))) +
  geom_bar(width = 0.72, fill = "#630084") +
  scale_x_discrete(name = "Optimal Si\u2013temperature lag (years)") +
  scale_y_continuous(name = "Number of sites",
                     expand = expansion(mult = c(0, 0.08))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey93"),
    plot.margin        = ggplot2::margin(0, 4, 0, 4)
  )

# =============================================================
# PANEL C — SMK tau by geology class + summer vs winter
# =============================================================
smk_geo         <- readRDS("data/smk/smk_geo.rds")
smk_seasonal_geo <- readRDS("data/smk/smk_seasonal_geo.rds")

# Full-period tau per site-geology (for ordering)
geo_order <- smk_geo %>%
  filter(variable == "Si_mgL", !is.na(tau), !is.na(geology_cat),
         geology_cat != "Other") %>%
  group_by(geology_cat) %>%
  summarise(med_tau = median(tau, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_tau) %>%
  pull(geology_cat) %>%
  as.character()

# Seasonal tau for violin
season_dat <- smk_seasonal_geo %>%
  filter(!is.na(tau), !is.na(geology_cat), geology_cat != "Other") %>%
  mutate(
    geology_cat = factor(as.character(geology_cat), levels = geo_order),
    season      = factor(season, levels = c("Summer (JJA)", "Winter (DJF)"))
  )

# Wilcoxon by geology class
wilcox_by_geo <- season_dat %>%
  select(site_no, geology_cat, season, tau) %>%
  pivot_wider(names_from = season, values_from = tau, values_fn = mean) %>%
  filter(!is.na(`Summer (JJA)`), !is.na(`Winter (DJF)`)) %>%
  group_by(geology_cat) %>%
  summarise(
    n     = n(),
    p_val = tryCatch(
      wilcox.test(`Summer (JJA)`, `Winter (DJF)`, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sig_label = case_when(p_val <= 0.001 ~ "***",
                          p_val <= 0.01  ~ "**",
                          p_val <= 0.05  ~ "*",
                          TRUE           ~ ""),
    ypos = 0.65   # annotation y position
  )

cat("\nWilcoxon results by geology:\n")
print(wilcox_by_geo)

# Overall Wilcoxon (all geologies combined) — for methods/results text
season_wide_all <- season_dat %>%
  select(site_no, geology_cat, season, tau) %>%
  pivot_wider(names_from = season, values_from = tau, values_fn = mean) %>%
  filter(!is.na(`Summer (JJA)`), !is.na(`Winter (DJF)`))

wilcox_overall <- wilcox.test(
  season_wide_all$`Summer (JJA)`,
  season_wide_all$`Winter (DJF)`,
  paired = TRUE, alternative = "greater"
)
cat("\nOverall Wilcoxon (summer > winter):\n")
cat("  W =", wilcox_overall$statistic,
    "  p =", format(wilcox_overall$p.value, digits = 4), "\n")

# Median tau per geology (for results text)
cat("\nMedian tau by geology (full period):\n")
smk_geo %>%
  filter(variable == "Si_mgL", !is.na(tau), !is.na(geology_cat),
         geology_cat != "Other") %>%
  group_by(geology_cat) %>%
  summarise(
    n       = n(),
    med_tau = round(median(tau, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(med_tau)) %>%
  print()

smk_geo <- readRDS("data/smk/smk_geo.rds")

geo_order <- smk_geo %>%
  filter(variable == "Si_mgL", !is.na(tau), 
         !is.na(geology_cat), geology_cat != "Other") %>%
  group_by(geology_cat) %>%
  summarise(med_tau = median(tau, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_tau) %>%
  pull(geology_cat) %>%
  as.character()

plot_dat <- smk_geo %>%
  filter(variable == "Si_mgL", !is.na(tau),
         !is.na(geology_cat), geology_cat != "Other") %>%
  mutate(geology_cat = factor(as.character(geology_cat), levels = geo_order))

geo_colors <- c(
  "Igneous"       = "#4A85B8",
  "Metamorphic"   = "#5A9E6F",
  "Sandstone"     = "#8B72B0",
  "Siltstone"     = "#C9B84C",
  "Clay/Mudstone" = "#C07B7B"
)

pC <- ggplot(plot_dat, aes(x = geology_cat, y = tau,
                           fill = geology_cat, color = geology_cat)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.45, color = "grey40") +
  geom_violin(trim = TRUE, scale = "width",
              alpha = 0.70, width = 0.68, color = NA) +
  geom_boxplot(width = 0.10, outlier.shape = NA,
               fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = geo_colors) +
  scale_color_manual(values = geo_colors, guide = "none") +
  scale_y_continuous(limits = c(-0.55, 0.70),
                     breaks = seq(-0.5, 0.5, by = 0.25),
                     name   = "SMK \u03c4 (WY2000\u2013WY2025)") +
  labs(x = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(color = "grey93"),
    panel.grid.minor = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.x      = element_text(size = 8.5, color = "grey20"),
    axis.title.y     = element_text(size = 9, margin = ggplot2::margin(r = 8)),
    plot.margin      = ggplot2::margin(0, 4, 4, 4)
  )

# =============================================================
# COMBINE
# =============================================================
fig2 <- (pA / pB / pC +
           plot_layout(heights = c(1.1, 0.9, 1.2))) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

ggsave(file.path(FIG_DIR, "figure2.pdf"), fig2, width = 7, height = 10)
ggsave(file.path(FIG_DIR, "figure2.png"), fig2, width = 7, height = 10, dpi = 300)
cat("Saved: figure2.pdf / figure2.png\n")

fig2
