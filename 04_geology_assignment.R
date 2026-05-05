  library(sf)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(readr)

GEOL_SHP <- "/Users/u1322101/Library/CloudStorage/OneDrive-UniversityofUtah/Documents_OneDrive/Research Projects/silica_temp/gagesII/geol_poly/geol_poly.shp"
WS_RDS   <- "/Users/u1322101/data/smk/ws_all_combined.rds"   
SMK_FILE <- "data/smk/smk_full.rds"
SEA_FILE <- "data/smk/smk_seasonal.rds"
OUT_DIR  <- "data/smk"

smk_sites <- readRDS(SMK_FILE) %>%
  distinct(site_no) %>%
  pull(site_no)


ws <- ws %>%
  mutate(GAGE_ID = str_pad(as.character(GAGE_ID), 8, pad = "0")) %>%
  filter(GAGE_ID %in% smk_sites) %>%
  st_transform(5070) %>%          
  select(GAGE_ID, geometry)


geol <- st_read(GEOL_SHP, quiet = TRUE)
stopifnot("LITH62" %in% names(geol))

geol <- geol %>%
  st_transform(5070) %>%
  st_make_valid() %>%
  select(LITH62, geometry) %>%
  mutate(
    LITH62 = as.character(LITH62),
    LITH62 = if_else(is.na(LITH62) | str_squish(LITH62) == "", "unknown", LITH62)
  )

classify_lith62 <- function(lith) {
  lith_low <- str_squish(str_to_lower(as.character(lith)))
  lith_low <- if_else(is.na(lith_low), "unknown", lith_low)

  case_when(
    str_detect(lith_low, paste(c(
      "limestone", "dolostone", "carbonate", "marble",
      "calcarenite", "travertine", "calc-silicate"
    ), collapse = "|")) ~ "Carbonate",

    str_detect(lith_low, paste(c(
      "metamorph", "gneiss", "schist", "quartzite", "serpentinite",
      "slate", "phyllite", "amphibolite", "migmatite", "greenstone",
      "granulite", "mylonite", "hornfels", "granofels", "phyllonite",
      "tectonite", "meta-argillite", "metasedimentary", "metavolcanic"
    ), collapse = "|")) ~ "Metamorphic",

    str_detect(lith_low, paste(c(
      "igneous", "granite", "granit", "granodiorite", "granitoid",
      "monzonite", "syenite", "tonalite", "gabbro", "diorite", "norite",
      "anorthosite", "diabase", "dolerite", "basalt", "rhyolite", "dacite",
      "andesite", "latite", "trachyte", "phonolite", "plutonic", "intrusive",
      "volcanic", "tuff", "ignimbrite", "lava", "ultramafic", "peridotite",
      "dunite", "trondhjemite", "pegmatite", "pyroclastic", "pyroxenite"
    ), collapse = "|")) ~ "Igneous",

    str_detect(lith_low, paste(c(
      "sandstone", "graywacke", "greywacke", "arenite", "arkose",
      "conglomerate", "coarse-grained mixed clastic"
    ), collapse = "|")) ~ "Sandstone",

    str_detect(lith_low, paste(c(
      "siltstone", "\\bsilt\\b", "fine-grained mixed clastic",
      "medium-grained mixed clastic", "clastic rock",
      "mixed volcanic/clastic"
    ), collapse = "|")) ~ "Siltstone",

    str_detect(lith_low, paste(c(
      "mudstone", "\\bclay\\b", "\\bmud\\b", "shale", "black shale",
      "argillite", "claystone", "chert", "novaculite", "iron formation"
    ), collapse = "|")) ~ "Clay/Mudstone",

    str_detect(lith_low, "\\bwater\\b|\\bice\\b") ~ NA_character_,

    TRUE ~ "Other"
  )
}

geol <- geol %>%
  mutate(geology_cat = classify_lith62(LITH62))

ws_geol <- st_intersection(ws, geol %>% select(LITH62, geology_cat)) %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry()


geol_frac <- ws_geol %>%
  filter(!is.na(geology_cat)) %>%       # drop water/ice
  group_by(GAGE_ID, geology_cat) %>%
  summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") %>%
  group_by(GAGE_ID) %>%
  mutate(
    total_area_m2 = sum(area_m2),
    frac_area     = area_m2 / total_area_m2
  ) %>%
  ungroup()


geo_dominant <- geol_frac %>%
  filter(geology_cat != "Carbonate") %>%
  arrange(GAGE_ID, desc(frac_area)) %>%
  group_by(GAGE_ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    site_no       = GAGE_ID,
    geology_cat   = geology_cat,
    frac_dominant = frac_area
  )

carbonate_frac <- geol_frac %>%
  filter(geology_cat == "Carbonate") %>%
  transmute(site_no = GAGE_ID, frac_carbonate = frac_area)

geol_wide <- geol_frac %>%
  select(GAGE_ID, geology_cat, frac_area) %>%
  pivot_wider(
    names_from   = geology_cat,
    values_from  = frac_area,
    values_fill  = 0,
    names_prefix = "frac_"
  ) %>%
  rename(site_no = GAGE_ID)


geo_site <- geo_dominant %>%
  left_join(carbonate_frac, by = "site_no") %>%
  left_join(geol_wide,      by = "site_no") %>%
  mutate(
    frac_carbonate = replace_na(frac_carbonate, 0),
    geology_cat = factor(
      geology_cat,
      levels = c("Igneous", "Metamorphic", "Sandstone",
                 "Siltstone", "Clay/Mudstone", "Other")
    )
  )

saveRDS(geo_site, file.path(OUT_DIR, "geo_site.rds"))

smk_full <- readRDS(SMK_FILE)

smk_geo <- smk_full %>%
  left_join(
    geo_site %>% select(site_no, geology_cat, frac_dominant, frac_carbonate),
    by = "site_no"
  )

saveRDS(smk_geo, file.path(OUT_DIR, "smk_geo.rds"))
smk_seasonal <- readRDS(SEA_FILE)
smk_seasonal_geo <- smk_seasonal %>%
  left_join(
    geo_site %>% select(site_no, geology_cat, frac_dominant, frac_carbonate),
    by = "site_no"
  )
saveRDS(smk_seasonal_geo, file.path(OUT_DIR, "smk_seasonal_geo.rds"))
smk_rolling <- readRDS(file.path(OUT_DIR, "smk_rolling.rds"))
smk_rolling_geo <- smk_rolling %>%
  left_join(
    geo_site %>% select(site_no, geology_cat, frac_dominant, frac_carbonate),
    by = "site_no"
  )
saveRDS(smk_rolling_geo, file.path(OUT_DIR, "smk_rolling_geo.rds"))

