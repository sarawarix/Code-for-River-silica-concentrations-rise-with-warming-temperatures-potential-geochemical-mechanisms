library(tidyverse)
MW_SiO2 <- 60.0843
MW_Si   <- 28.0855
OUT_DIR <- "data/wqp_download"

BASIN_FILE <- "/Users/u1322101/Library/CloudStorage/OneDrive-UniversityofUtah/Documents_OneDrive/Research Projects/silica_temp/gagesII/Dataset1_BasinID/BasinID.txt"

wqp_clean <- wqp_screened %>%
  select(Location_Identifier, Activity_StartDate,
         Result_Characteristic, Result_Measure, Result_MeasureUnit) %>%
  mutate(
    Date           = as.Date(Activity_StartDate),
    site_no        = str_remove(Location_Identifier, "^USGS-"),
    Result_Measure = suppressWarnings(as.numeric(Result_Measure)),
    char_lower     = tolower(Result_Characteristic)
  ) %>%
  filter(!is.na(Result_Measure), Result_Measure >= 0, !is.na(Date)) %>%
  mutate(
    value_mgL = case_when(
      char_lower == "silica"  & Result_MeasureUnit == "mg/L"           ~ Result_Measure,
      char_lower == "silicon" & Result_MeasureUnit == "ug/L"           ~
        Result_Measure / 1000 * (MW_SiO2 / MW_Si),
      char_lower == "ph"      & Result_MeasureUnit == "standard units" ~ Result_Measure,
      TRUE ~ NA_real_
    ),
    param = case_when(
      char_lower %in% c("silica", "silicon") ~ "Si_mgL",
      char_lower == "ph"                     ~ "pH"
    )
  ) %>%
  filter(!is.na(value_mgL), !is.na(param))

wqp_wide <- wqp_clean %>%
  group_by(site_no, Date, param) %>%
  summarise(value = median(value_mgL, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = value)

cq <- wqp_wide %>%
  inner_join(nwis_slim, by = c("site_no", "Date")) %>%
  filter(
    !is.na(Flow),
    Date >= as.Date("2000-01-01"),
    Date <= as.Date("2025-09-30")
  )



basin <- read.csv(BASIN_FILE, stringsAsFactors = FALSE, strip.white = TRUE) %>%
  mutate(STAID = str_pad(as.character(STAID), 8, pad = "0")) %>%
  select(STAID, HUC02, STATE, CLASS, AGGECOREGION)

cq <- cq %>%
  mutate(
    Q_Ls        = Flow * 28.3168,
    Si_flux_mgs = Si_mgL * Q_Ls
  ) %>%
  left_join(basin, by = c("site_no" = "STAID")) %>%
  filter(Flow > 0, Si_mgL > 0, Si_mgL <= 50) %>%
  select(site_no, Date, Si_mgL, Ca_mgL, Na_mgL, pH,
         Flow, Q_Ls, Wtemp, Si_flux_mgs,
         HUC02, STATE, CLASS, AGGECOREGION) %>%
  arrange(site_no, Date)

saveRDS(cq, file.path(OUT_DIR, "cq_master.rds"))
write_csv(cq, file.path(OUT_DIR, "cq_master.csv"))
