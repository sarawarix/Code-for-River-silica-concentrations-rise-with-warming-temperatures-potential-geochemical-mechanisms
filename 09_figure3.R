library(tidyverse)
library(stringr)
library(readr)
library(patchwork)

base_dir <- enc2utf8("/Users/u1322101/Documents/v2_manuscript_models")
scenarios <- c("run1_T_cT","run2_T_dT","run3_K_cT","run4_K_dT")

R <- 8.314462618
MW_SiO2 <- 60.0843  # g/mol

T_pts <- c(0, 25)
logK_alb  <- c( 3.2730,  2.7645)
logK_kaol <- c( -9.0182, -6.8101)

logK_comb_pts <- 2*logK_alb + logK_kaol   # old way
log10K_comb <- function(T_C) approx(T_pts, logK_comb_pts, xout = T_C, rule = 2)$y

list_ordered <- function(dirpath, pattern) {
  files <- list.files(dirpath, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop(paste("No files found for pattern:", pattern, "in", dirpath))
  files[order(as.numeric(str_extract(basename(files), "\\d+")))]
}

read_x_from_last_totconread_x_from_last_totcon <- function(run_path) {
  f <- tail(list_ordered(run_path, "^totcon\\d+\\.out$"), 1)
  
  ln <- readLines(f, warn = FALSE, encoding = "latin1")
  ln <- iconv(ln, from = "latin1", to = "UTF-8", sub = "")
  
  header_idx <- which(grepl("^\\s*Distance\\b", ln))[1]
  header_tokens <- strsplit(trimws(ln[header_idx]), "\\s+")[[1]]
  header_tokens[1] <- "x_node_m"
  
  df <- read_table(
    f, skip = header_idx, comment = "#",
    col_names = header_tokens, col_types = cols(.default = col_double()),
    progress = FALSE
  )
  
  nm <- names(df)
  if ("SiO2(aq)" %in% nm) df <- dplyr::rename(df, SiO2_tot = `SiO2(aq)`)
  
  df %>%
    arrange(x_node_m) %>%
    mutate(node_index = row_number()) %>%
    select(node_index, x_node_m, any_of("SiO2_tot"))
}
read_from_speciation_file <- function(spec_file) {
  ln <- readLines(spec_file, warn = FALSE)
  loc_idx <- grep("GRID LOCATION", ln)
  if (!length(loc_idx)) return(tibble())
  
  block_starts <- loc_idx
  block_ends   <- c(loc_idx[-1] - 1, length(ln))
  
  parse_species_row <- function(table_lines, spec_pat) {
    i <- grep(paste0("^\\s*", spec_pat, "\\s"), table_lines)
    if (!length(i)) return(tibble(m = NA_real_, loga = NA_real_))
    tokens <- strsplit(trimws(table_lines[i[1]]), "\\s+")[[1]]
    tibble(
      m    = suppressWarnings(as.numeric(tokens[2])),
      loga = suppressWarnings(as.numeric(tokens[3]))
    )
  }
  
  parse_SI <- function(block, mineral_name) {
    i <- grep(paste0("^\\s*", mineral_name, "\\b"), block)
    if (!length(i)) return(NA_real_)
    nums <- str_extract_all(block[i[1]], "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?")[[1]]
    if (!length(nums)) return(NA_real_)
    suppressWarnings(as.numeric(nums[1]))
  }
  
  purrr::map2_dfr(block_starts, block_ends, function(s, e) {
    block <- ln[s:e]
    node <- as.integer(str_match(block[1], "GRID LOCATION:\\s*(\\d+)")[,2])
    
    t_idx <- grep("Temperature \\(C\\)", block)
    T_C <- if (length(t_idx)) as.numeric(str_extract(block[t_idx[1]], "[-0-9.]+$")) else NA_real_
    
    hdr_idx <- grep(
      "^\\s*Species\\s+Molality\\s+Activity\\s+Molality\\s+Activity\\s+Coefficient\\s+Type",
      block
    )[1]
    
    out <- tibble(
      node_index = node,
      T_C        = T_C,
      Na_m       = NA_real_, loga_Na   = NA_real_,
      H_m        = NA_real_, loga_H    = NA_real_,
      SiO2_m     = NA_real_, loga_SiO2 = NA_real_,
      SI_kaol    = parse_SI(block, "Kaolinite"),
      SI_alb     = parse_SI(block, "Albite")
    )
    
    if (is.na(hdr_idx)) return(out)
    
    table_lines <- block[(hdr_idx + 1):length(block)]
    blank_idx <- which(trimws(table_lines) == "")
    if (length(blank_idx)) table_lines <- table_lines[1:(blank_idx[1]-1)]
    
    Na <- parse_species_row(table_lines, "Na\\+")
    H  <- parse_species_row(table_lines, "H\\+")
    Si <- parse_species_row(table_lines, "SiO2\\(aq\\)")
    
    out %>%
      mutate(
        Na_m      = Na$m,  loga_Na   = Na$loga,
        H_m       = H$m,   loga_H    = H$loga,
        SiO2_m    = Si$m,  loga_SiO2 = Si$loga
      )
  }) %>% arrange(node_index)
}

compute_profiles_fast <- function(run_path) {
  x_df <- read_x_from_last_totcon(run_path)
  
  spec_file <- tail(list_ordered(run_path, "^speciation\\d+\\.out$"), 1)
  a_df <- read_from_speciation_file(spec_file)
  
  x_df %>%
    left_join(a_df, by = "node_index") %>%
    mutate(
      T_K = T_C + 273.15,
      
      log10Q = 2*loga_Na + 4*loga_SiO2 - 2*loga_H,
      log10K = log10K_comb(T_C),
      
      dG_Jmol  = 2.303 * R * T_K * (log10Q - log10K),
      dG_kJmol = dG_Jmol / 1000,
      
      SiO2_molL = dplyr::coalesce(SiO2_tot, SiO2_m),
      SiO2_mgL  = SiO2_molL * MW_SiO2 * 1000
    )
}

profiles <- purrr::map_dfr(scenarios, function(scn) {
  run_path <- file.path(base_dir, scn)
  compute_profiles_fast(run_path) %>% mutate(scenario = scn)
}) %>%
  mutate(
    # UPDATED to match new names
    process = case_when(
      scenario %in% c("run1_T_cT", "run2_T_dT") ~ "Thermodynamic Eq.",
      scenario %in% c("run3_K_cT", "run4_K_dT") ~ "Kinetic Buffering"
    ),
    temp = case_when(
      scenario %in% c("run2_T_dT", "run4_K_dT") ~ "Warming",
      scenario %in% c("run1_T_cT", "run3_K_cT") ~ "Constant temperature"
    )
  )

profiles_solid  <- profiles %>% filter(temp == "Warming")
profiles_dashed <- profiles %>% filter(temp == "Constant temperature")

# -----------------------------
# Panel A: [SiO2(aq)]
# -----------------------------
pA <- ggplot() +
  geom_line(
    data = profiles_solid,
    aes(x = x_node_m, y = SiO2_mgL, color = process, group = interaction(process, scenario)),
    linewidth = 1,
    linetype = "solid"
  ) +
  geom_line(
    data = profiles_dashed,
    aes(x = x_node_m, y = SiO2_mgL, color = process, group = interaction(process, scenario)),
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_color_manual(values = c("Thermodynamic Eq." = "grey60",
                                "Kinetic Buffering" = "black")) +
  scale_y_continuous(
    name = expression("[SiO"[2]*"(aq)] (mg/L)"),
    sec.axis = sec_axis(~ . / (MW_SiO2 * 1000),
                        name = expression("[SiO"[2]*"(aq)] (mol/L)"))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# -----------------------------
# Panel B: Saturation Index (same scale; correct grouping; same linetype workaround)
# -----------------------------
# clearer differentiation:
# - color = mineral (2 colors)
# - linetype = scenario (4 linetypes)
# - optional: linewidth by process (T vs K) for an extra cue

si_long <- profiles %>%
  select(x_node_m, process, scenario, SI_kaol, SI_alb) %>%
  pivot_longer(c(SI_kaol, SI_alb), names_to = "mineral", values_to = "SI") %>%
  mutate(
    mineral  = recode(mineral, SI_kaol = "Kaolinite", SI_alb = "Albite"),
    scenario = factor(scenario, levels = scenarios),
    process  = factor(process, levels = c("Thermodynamic Eq.", "Kinetic Buffering"))
  )

pB <- ggplot(
  si_long,
  aes(
    x = x_node_m,
    y = SI,
    color = mineral,
    linetype = scenario,
    group = interaction(mineral, scenario),
    linewidth = process
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line() +
  scale_linewidth_manual(values = c("Thermodynamic Eq." = 0.8, "Kinetic Buffering" = 1.2)) +
  # if you want specific linetypes per run (optional)
  # scale_linetype_manual(values = c("run1_T_cT"="solid","run2_T_dT"="dashed","run3_K_cT"="dotdash","run4_K_dT"="twodash")) +
  labs(x = NULL, y = "Saturation index") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")


# -----------------------------
# Panel C: ΔG (same linetype workaround)
# -----------------------------
pC <- ggplot() +
  geom_line(
    data = profiles_solid,
    aes(x = x_node_m, y = dG_kJmol, color = process, group = interaction(process, scenario)),
    linewidth = 1,
    linetype = "solid"
  ) +
  geom_line(
    data = profiles_dashed,
    aes(x = x_node_m, y = dG_kJmol, color = process, group = interaction(process, scenario)),
    linewidth = 1,
    linetype = "dashed"
  ) +
  #coord_cartesian(ylim = c(-117.5, -110)) +
  scale_color_manual(values = c("Thermodynamic Eq." = "grey60",
                                "Kinetic Buffering" = "black")) +
  labs(x = "Distance (m)", y = expression(Delta*G[comb]~"(kJ/mol)")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# -----------------------------
# Combine
# -----------------------------
(pA / pB) +
  plot_layout(heights = c(1, 1.2, 1))


#### stats for manuscript
#thermodynamic scenarios
thermo <- profiles %>% filter(process == "Thermodynamic Eq.")
thermo_summary <- thermo %>%
  group_by(scenario, temp) %>%
  summarise(
    SiO2_min = min(SiO2_mgL, na.rm = TRUE),
    SiO2_max = max(SiO2_mgL, na.rm = TRUE),
    SiO2_increase = max(SiO2_mgL, na.rm = TRUE) - min(SiO2_mgL, na.rm = TRUE),
    .groups = "drop"
  )
increase_warming  <- thermo_summary %>% filter(temp == "Warming") %>% pull(SiO2_increase)
increase_constant <- thermo_summary %>% filter(temp == "Constant temperature") %>% pull(SiO2_increase)
pct_higher <- (increase_warming - increase_constant) / increase_constant * 100

print(thermo_summary)
cat("Warming increase is:", round(pct_higher, 1), "% higher\n")

#thermodynamic saturation
thermo_SI_summary <- si_long %>%
  filter(process == "Thermodynamic Eq.") %>%
  group_by(mineral, scenario) %>%
  summarise(
    SI_start = first(SI),
    SI_end   = last(SI),
    SI_min   = min(SI, na.rm = TRUE),
    SI_max   = max(SI, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(temp = case_when(
    scenario == "run1_T_cT" ~ "Constant temperature",
    scenario == "run2_T_dT" ~ "Warming"
  ))

print(thermo_SI_summary)

#kinetic model
kinetic <- profiles %>% filter(process == "Kinetic Buffering")

kinetic_summary <- kinetic %>%
  group_by(scenario, temp) %>%
  summarise(
    SiO2_start  = first(SiO2_mgL),
    SiO2_end    = last(SiO2_mgL),
    SiO2_max    = max(SiO2_mgL, na.rm = TRUE),
    SiO2_min    = min(SiO2_mgL, na.rm = TRUE),
    SiO2_increase = max(SiO2_mgL, na.rm = TRUE) - min(SiO2_mgL, na.rm = TRUE),
    x_at_max    = x_node_m[which.max(SiO2_mgL)],  # distance where plateau begins
    .groups = "drop"
  )

kinetic_SI_summary <- si_long %>%
  filter(process == "Kinetic Buffering") %>%
  group_by(mineral, scenario) %>%
  summarise(
    SI_start = first(SI),
    SI_end   = last(SI),
    SI_min   = min(SI, na.rm = TRUE),
    SI_max   = max(SI, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(temp = case_when(
    scenario == "run3_K_cT" ~ "Constant temperature",
    scenario == "run4_K_dT" ~ "Warming"
  ))

print(kinetic_summary)
print(kinetic_SI_summary)

const_kinetic <- profiles %>% 
  filter(scenario == "run3_K_cT") %>%
  arrange(x_node_m) %>%
  mutate(
    SiO2_change = abs(SiO2_mgL - lag(SiO2_mgL)),  # change between nodes
  )

const_kinetic %>% 
  select(x_node_m, SiO2_mgL, SiO2_change) %>% 
  print(n = Inf)

