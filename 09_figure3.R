library(tidyverse)
library(stringr)
library(readr)
library(patchwork)

R_gas   <- 8.314462618
MW_SiO2 <- 60.0843

base_dir  <- enc2utf8("/Users/u1322101/Library/CloudStorage/OneDrive-UniversityofUtah/Documents_OneDrive/Research Projects/silica_temp/Silica_Temp_ERL/codeDump/v4_manuscript_models")
scenarios <- c("run1_T_cT", "run2_T_dT_small", "run3_T_dT_big",
               "run4_K_cT", "run5_K_dT_small", "run6_K_dT_big")

T_pts        <- c(0, 25)
logK_alb     <- c( 3.2730,  2.7645)
logK_kaol    <- c(-9.0182, -6.8101)
logK_comb_pts <- 2 * logK_alb + logK_kaol
log10K_comb  <- function(T_C) approx(T_pts, logK_comb_pts, xout = T_C, rule = 2)$y

list_ordered <- function(dirpath, pattern) {
  files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop(paste("No files found:", pattern, "in", dirpath))
  files[order(as.numeric(str_extract(basename(files), "\\d+")))]
}

read_x_from_last_totcon <- function(run_path) {
  f  <- tail(list_ordered(run_path, "^totcon\\d+\\.out$"), 1)
  ln <- readLines(f, warn = FALSE, encoding = "latin1")
  ln <- iconv(ln, from = "latin1", to = "UTF-8", sub = "")
  
  header_idx    <- which(grepl("^\\s*Distance\\b", ln))[1]
  header_tokens <- strsplit(trimws(ln[header_idx]), "\\s+")[[1]]
  header_tokens[1] <- "x_node_m"
  
  df <- read_table(
    f, skip = header_idx, comment = "#",
    col_names = header_tokens, col_types = cols(.default = col_double()),
    progress  = FALSE
  )
  
  if ("SiO2(aq)" %in% names(df)) df <- dplyr::rename(df, SiO2_tot = `SiO2(aq)`)
  
  df %>%
    arrange(x_node_m) %>%
    mutate(node_index = row_number()) %>%
    select(node_index, x_node_m, any_of("SiO2_tot"))
}

read_from_speciation_file <- function(spec_file) {
  ln      <- readLines(spec_file, warn = FALSE)
  loc_idx <- grep("GRID LOCATION", ln)
  if (!length(loc_idx)) return(tibble())
  
  parse_species_row <- function(table_lines, spec_pat) {
    i <- grep(paste0("^\\s*", spec_pat, "\\s"), table_lines)
    if (!length(i)) return(tibble(m = NA_real_, loga = NA_real_))
    tokens <- strsplit(trimws(table_lines[i[1]]), "\\s+")[[1]]
    tibble(m = suppressWarnings(as.numeric(tokens[2])),
           loga = suppressWarnings(as.numeric(tokens[3])))
  }
  
  parse_SI <- function(block, mineral_name) {
    i <- grep(paste0("^\\s*", mineral_name, "\\b"), block)
    if (!length(i)) return(NA_real_)
    nums <- str_extract_all(block[i[1]], "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?")[[1]]
    if (!length(nums)) return(NA_real_)
    suppressWarnings(as.numeric(nums[1]))
  }
  
  purrr::map2_dfr(loc_idx, c(loc_idx[-1] - 1, length(ln)), function(s, e) {
    block  <- ln[s:e]
    node   <- as.integer(str_match(block[1], "GRID LOCATION:\\s*(\\d+)")[, 2])
    t_idx  <- grep("Temperature \\(C\\)", block)
    T_C    <- if (length(t_idx)) as.numeric(str_extract(block[t_idx[1]], "[-0-9.]+$")) else NA_real_
    
    hdr_idx <- grep(
      "^\\s*Species\\s+Molality\\s+Activity\\s+Molality\\s+Activity\\s+Coefficient\\s+Type",
      block
    )[1]
    
    out <- tibble(
      node_index = node, T_C = T_C,
      Na_m = NA_real_, loga_Na   = NA_real_,
      H_m  = NA_real_, loga_H    = NA_real_,
      SiO2_m = NA_real_, loga_SiO2 = NA_real_,
      SI_kaol = parse_SI(block, "Kaolinite"),
      SI_alb  = parse_SI(block, "Albite")
    )
    
    if (is.na(hdr_idx)) return(out)
    
    table_lines <- block[(hdr_idx + 1):length(block)]
    blank_idx   <- which(trimws(table_lines) == "")
    if (length(blank_idx)) table_lines <- table_lines[1:(blank_idx[1] - 1)]
    
    Na <- parse_species_row(table_lines, "Na\\+")
    H  <- parse_species_row(table_lines, "H\\+")
    Si <- parse_species_row(table_lines, "SiO2\\(aq\\)")
    
    out %>% mutate(
      Na_m = Na$m, loga_Na   = Na$loga,
      H_m  = H$m,  loga_H    = H$loga,
      SiO2_m = Si$m, loga_SiO2 = Si$loga
    )
  }) %>% arrange(node_index)
}

compute_profiles_fast <- function(run_path) {
  x_df      <- read_x_from_last_totcon(run_path)
  spec_file <- tail(list_ordered(run_path, "^speciation\\d+\\.out$"), 1)
  a_df      <- read_from_speciation_file(spec_file)
  
  x_df %>%
    left_join(a_df, by = "node_index") %>%
    mutate(
      T_K       = T_C + 273.15,
      log10Q    = 2 * loga_Na + 4 * loga_SiO2 - 2 * loga_H,
      log10K    = log10K_comb(T_C),
      dG_Jmol   = 2.303 * R_gas * T_K * (log10Q - log10K),
      dG_kJmol  = dG_Jmol / 1000,
      SiO2_molL = dplyr::coalesce(SiO2_tot, SiO2_m),
      SiO2_mgL  = SiO2_molL * MW_SiO2 * 1000
    )
}

# ── Build profiles ────────────────────────────────────────────────
profiles <- purrr::map_dfr(scenarios, function(scn) {
  compute_profiles_fast(file.path(base_dir, scn)) %>% mutate(scenario = scn)
}) %>%
  mutate(
    process = case_when(
      scenario %in% c("run1_T_cT", "run2_T_dT_small", "run3_T_dT_big") ~ "Thermodynamic Eq.",
      scenario %in% c("run4_K_cT", "run5_K_dT_small", "run6_K_dT_big") ~ "Kinetic Buffering"
    ),
    temp = factor(case_when(
      scenario %in% c("run3_T_dT_big",   "run6_K_dT_big")   ~ "High warming",
      scenario %in% c("run2_T_dT_small", "run5_K_dT_small") ~ "Low warming",
      scenario %in% c("run1_T_cT",       "run4_K_cT")       ~ "Constant temperature"
    ), levels = c("High warming", "Low warming", "Constant temperature"))
  )

base_font <- "sans"

theme_manuscript <- function() {
  theme_bw(base_size = 10) +
    theme(
      text              = element_text(family = base_font),
      panel.grid.major  = element_line(color = "grey93"),
      panel.grid.minor  = element_blank(),
      axis.ticks        = element_line(color = "grey40", linewidth = 0.3),
      axis.text         = element_text(size = 8.5, color = "grey20"),
      axis.title        = element_text(size = 9),
      axis.title.y.right = element_text(color = "grey40", size = 9),
      axis.text.y.right  = element_text(color = "grey40"),
      axis.ticks.y.right = element_line(color = "grey40"),
      plot.margin       = ggplot2::margin(4, 4, 4, 4)
    )
}

lt_vals   <- c("High warming" = "solid", "Low warming" = "dashed", "Constant temperature" = "dotted")
lt_labels <- c("High warming (1.5°C/decade)", "Low warming (0.15°C/decade)", "Constant temp.")

annot_B <- profiles %>%
  filter(temp %in% c("High warming", "Constant temperature")) %>%
  group_by(process, temp) %>%
  summarise(y = max(SiO2_mgL, na.rm = TRUE), .groups = "drop")

pB <- ggplot(
  profiles,
  aes(
    x        = x_node_m,
    y        = SiO2_mgL,
    color    = process,
    linetype = temp,
    group    = interaction(process, temp)
  )
) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = c("Thermodynamic Eq." = "grey60", "Kinetic Buffering" = "black"),
    name   = NULL
  ) +
  scale_linetype_manual(
    values = lt_vals,
    labels = lt_labels,
    name   = NULL
  ) +
  scale_y_continuous(
    name = expression("[SiO"[2]*"(aq)]  (mg/L)"),
    sec.axis = sec_axis(
      ~ . / (MW_SiO2 * 1000),
      name   = expression("[SiO"[2]*"(aq)]  (mol/L)"),
      labels = scales::scientific
    )
  ) +
  annotate("text", x = 280,
           y = profiles$SiO2_mgL[profiles$scenario == "run6_K_dT_big"] %>% max(na.rm = TRUE),
           label = "High warming", hjust = 0, vjust = -0.4,
           size = 3, family = base_font) +
  annotate("text", x = 280,
           y = profiles$SiO2_mgL[profiles$scenario == "run4_K_cT"] %>% max(na.rm = TRUE),
           label = "Constant temp.", hjust = 0, vjust = -0.4,
           size = 3, family = base_font) +
  # Process labels at x = 10
  annotate("text", x = 10,
           y = mean(profiles$SiO2_mgL[profiles$process == "Kinetic Buffering"], na.rm = TRUE),
           label = "Kinetic", hjust = 0, vjust = -0.5,
           size = 3, family = base_font, fontface = "bold") +
  annotate("text", x = 10,
           y = mean(profiles$SiO2_mgL[profiles$process == "Thermodynamic Eq."], na.rm = TRUE),
           label = "Equilibrium", hjust = 0, vjust = 1.2,
           size = 3, family = base_font, color = "grey50") +
  labs(x = NULL) +
  theme_manuscript() +
  theme(legend.position = "none")

si_plot <- profiles %>%
  select(x_node_m, process, temp, scenario, SI_kaol, SI_alb) %>%
  pivot_longer(c(SI_kaol, SI_alb), names_to = "mineral", values_to = "SI") %>%
  mutate(
    mineral     = recode(mineral, SI_kaol = "Kaolinite", SI_alb = "Albite"),
    model       = if_else(process == "Kinetic Buffering", "Kinetic", "Equilibrium"),
    color_group = paste(mineral, model)
  )

si_colors <- c(
  "Kaolinite Kinetic"     = "#2266BB",
  "Kaolinite Equilibrium" = "#88CCFF",
  "Albite Kinetic"        = "#CC3333",
  "Albite Equilibrium"    = "#FFAAAA"
)

pC <- ggplot(
  si_plot,
  aes(
    x         = x_node_m,
    y         = SI,
    color     = color_group,
    linetype  = temp,
    linewidth = model,
    group     = interaction(color_group, temp)
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_line() +
  scale_color_manual(values = si_colors, name = NULL) +
  scale_linetype_manual(
    values = lt_vals,
    labels = lt_labels,
    name   = NULL
  ) +
  scale_linewidth_manual(
    values = c("Kinetic" = 1.0, "Equilibrium" = 0.7),
    guide  = "none"
  ) +
  annotate("text", x = 310, y =  7.4,  label = "Kaolinite Kinetic",
           color = "#2266BB", hjust = 1, size = 2.8, family = base_font) +
  annotate("text", x = 310, y =  0.2,  label = "Kaolinite Equilibrium",
           color = "#88CCFF", hjust = 1, size = 2.8, family = base_font) +
  annotate("text", x = 310, y = -0.2,  label = "Albite Kinetic",
           color = "#CC3333", hjust = 1, size = 2.8, family = base_font) +
  annotate("text", x = 310, y = -4.8,  label = "Albite Equilibrium",
           color = "#FFAAAA", hjust = 1, size = 2.8, family = base_font) +
  labs(x = "Distance (m)", y = "Saturation index") +
  guides(
    color    = guide_legend(override.aes = list(linewidth = 0.8)),
    linetype = guide_legend(
      override.aes = list(linewidth = 0.8, color = "grey30"),
      title = NULL
    )
  ) +
  theme_manuscript() +
  theme(
    legend.position  = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.text      = element_text(size = 8)
  )

pB / pC +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12, family = base_font)
    )
  )
