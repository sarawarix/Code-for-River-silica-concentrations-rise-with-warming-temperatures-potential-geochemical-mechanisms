library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)

FIG_DIR <- "data/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ts_dat  <- readRDS("data/smk/ts_dat.rds")
ts_plot <- ts_dat %>% filter(wy >= 2000, wy <= 2024)

si_lm   <- lm(med_si   ~ wy, data = ts_plot)
tair_lm <- lm(med_tair ~ wy, data = ts_plot)

si_slope <- coef(si_lm)["wy"]
si_p     <- summary(si_lm)$coefficients["wy", "Pr(>|t|)"]
t_slope  <- coef(tair_lm)["wy"]
t_p      <- summary(tair_lm)$coefficients["wy", "Pr(>|t|)"]

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

si_label   <- sprintf("Si: %+.3f mg/L yr\u207b\u00b9,  p = %.4f", si_slope, si_p)
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
  geom_smooth(aes(y = tair_scaled, color = "Air temperature"),
              method = "lm", se = FALSE, linewidth = 0.7, linetype = "dashed") +
  geom_smooth(aes(y = med_si, color = "Silica"),
              method = "lm", se = FALSE, linewidth = 0.7, linetype = "dashed") +
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

optimal_lags <- readRDS("data/smk/optimal_lags.rds")

pos_dat <- optimal_lags %>%
  filter(direction == "positive") %>%
  mutate(
    sig_label = if_else(sig, "p \u2264 0.05", "p > 0.05"),
    sig_label = factor(sig_label, levels = c("p > 0.05", "p \u2264 0.05"))
  )

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
smk_rolling_geo     <- readRDS("data/smk/smk_rolling_geo.rds")
smk_seasonal_geo    <- readRDS("data/smk/smk_seasonal_geo.rds")

plot_dat <- smk_rolling_geo %>%
  filter(variable == "Si_mgL", wy_start == 2013, wy_end == 2022,
         !is.na(tau), !is.na(geology_cat), geology_cat != "Other") %>%
  mutate(
    lith_binary = case_when(
      geology_cat %in% c("Igneous", "Metamorphic")                  ~ "Crystalline",
      geology_cat %in% c("Sandstone", "Siltstone", "Clay/Mudstone") ~ "Sedimentary",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(lith_binary))

wt_geo <- wilcox.test(tau ~ lith_binary, data = plot_dat, alternative = "greater")

pC_geo <- ggplot(plot_dat, aes(x = lith_binary, y = tau, fill = lith_binary)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.45, color = "grey40") +
  geom_violin(trim = TRUE, scale = "width",
              alpha = 0.70, width = 0.68, color = NA) +
  geom_boxplot(width = 0.10, outlier.shape = NA,
               fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Crystalline" = "#8B72B0",
                               "Sedimentary" = "#C2AF4C")) +
  scale_y_continuous(limits = c(-0.55, 0.70),
                     breaks = seq(-0.5, 0.5, by = 0.25),
                     name   = "SMK \u03c4 (WY2013\u20132022)") +
  annotate("text", x = 1.5, y = 0.65,
           label = sprintf("Wilcoxon p = %.3f", wt_geo$p.value),
           size = 3.2, color = "grey20") +
  labs(x = NULL, fill = NULL, title = "Lithology") +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(color = "grey93"),
    panel.grid.minor = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.x      = element_text(size = 9, color = "grey20"),
    axis.title.y     = element_text(size = 9, margin = ggplot2::margin(r = 8)),
    plot.title       = element_text(size = 9, hjust = 0.5, color = "grey30"),
    plot.margin      = ggplot2::margin(0, 2, 4, 4)
  )

season_dat <- smk_seasonal_geo %>%
  filter(!is.na(tau), !is.na(geology_cat), geology_cat != "Other",
         season %in% c("Summer (JJA)", "Winter (DJF)")) %>%
  mutate(season = factor(season, levels = c("Summer (JJA)", "Winter (DJF)")))

season_wide <- season_dat %>%
  select(site_no, season, tau) %>%
  pivot_wider(names_from = season, values_from = tau, values_fn = mean) %>%
  filter(!is.na(`Summer (JJA)`), !is.na(`Winter (DJF)`))

wt_season <- wilcox.test(
  season_wide$`Summer (JJA)`,
  season_wide$`Winter (DJF)`,
  paired = TRUE, alternative = "greater"
)

pC_season <- ggplot(season_dat, aes(x = season, y = tau, fill = season)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.45, color = "grey40") +
  geom_violin(trim = TRUE, scale = "width",
              alpha = 0.70, width = 0.68, color = NA) +
  geom_boxplot(width = 0.10, outlier.shape = NA,
               fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Summer (JJA)" = "#E07B39",
                               "Winter (DJF)"  = "#4A85B8")) +
  scale_y_continuous(limits = c(-0.55, 0.70),
                     breaks = seq(-0.5, 0.5, by = 0.25),
                     name   = NULL) +
  annotate("text", x = 1.5, y = 0.65,
           label = sprintf("Wilcoxon p < 0.001"),
           size = 3.2, color = "grey20") +
  labs(x = NULL, fill = NULL, title = "Season") +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(color = "grey93"),
    panel.grid.minor = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.x      = element_text(size = 9, color = "grey20"),
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    plot.title       = element_text(size = 9, hjust = 0.5, color = "grey30"),
    plot.margin      = ggplot2::margin(0, 4, 4, 2)
  )

pC <- pC_geo | pC_season

fig2 <- (pA / pB / pC +
           plot_layout(heights = c(1.1, 0.9, 1.2))) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12))
  )

fig2
