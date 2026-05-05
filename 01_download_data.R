library(tidyverse)
library(dataRetrieval)

#script to download silica data for gages sites

MIN_SI_SAMPLES <- 30
START_DATE     <- "2000-01-01"   
END_DATE       <- "2025-09-30"             
OUT_DIR        <- "data/wqp_download"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BASIN_FILE <- "/Users/u1322101/Library/CloudStorage/OneDrive-UniversityofUtah/Documents_OneDrive/Research Projects/silica_temp/gagesII/Dataset1_BasinID/BasinID.txt"

basin <- read.csv(BASIN_FILE, stringsAsFactors = FALSE, strip.white = TRUE) %>%
  mutate(
    STAID  = str_pad(as.character(STAID), 8, pad = "0"),
    WQP_ID = paste0("USGS-", STAID)
  )

cat("Total GAGES-II sites:", nrow(basin), "\n")
all_wqp_ids <- basin$WQP_ID

wqp_chars <- list(
  Si = c("Silica", "Silicon"),
  pH = c("pH"))

all_chars <- unlist(wqp_chars, use.names = FALSE)

coerce_to_character <- function(df) {
  df %>% mutate(across(everything(), as.character))
}


batch_size   <- 100
site_batches <- split(all_wqp_ids, ceiling(seq_along(all_wqp_ids) / batch_size))

si_counts <- map_dfr(seq_along(site_batches), function(i) {
  batch <- site_batches[[i]]
  if (i %% 10 == 0) cat("  Screening batch", i, "of", length(site_batches), "\n")

  result <- tryCatch(
    whatWQPdata(
      siteNumber         = batch,
      characteristicName = c("Silica", "Silicon"),
      startDateLo        = START_DATE
    ),
    error = function(e) {
      cat("  Batch", i, "ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(result) || nrow(result) == 0) return(NULL)
  result
})


eligible_sites <- si_counts %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarise(n_si = sum(resultCount, na.rm = TRUE), .groups = "drop") %>%
  filter(n_si >= MIN_SI_SAMPLES) %>%
  pull(MonitoringLocationIdentifier)

saveRDS(eligible_sites, file.path(OUT_DIR, "eligible_sites.rds"))



wqp_dir <- file.path(OUT_DIR, "wqp_by_site")
dir.create(wqp_dir, showWarnings = FALSE)

for (site in eligible_sites) {
  out_file <- file.path(wqp_dir, paste0(str_remove(site, "^USGS-"), ".rds"))
  if (file.exists(out_file)) next  


  site_data <- map_dfr(all_chars, function(char) {
    result <- tryCatch(
      readWQPqw(
        siteNumbers = site,
        parameterCd = char,
        startDate   = START_DATE,
        endDate     = END_DATE,
        legacy      = FALSE     
      ),
      error = function(e) {
        cat("    ", char, "- ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(result) || nrow(result) == 0) return(NULL)
    cat("    ", char, "- n =", nrow(result), "\n")
    coerce_to_character(result)
  })

  if (is.null(site_data) || nrow(site_data) == 0) {
    cat("    no data\n")
    saveRDS(NULL, out_file)
    next
  }

  saveRDS(site_data, out_file)
}

site_files <- list.files(wqp_dir, pattern = "\\.rds$", full.names = TRUE)
has_data   <- map_lgl(site_files, ~ !is.null(readRDS(.x)))

wqp_raw <- map_dfr(site_files[has_data], readRDS)
saveRDS(wqp_raw, file.path(OUT_DIR, "wqp_raw_all.rds"))

wqp_screened <- wqp_raw %>%
  filter(Location_Identifier %in% eligible_sites)

saveRDS(wqp_screened, file.path(OUT_DIR, "wqp_screened_all.rds"))

usgs_sites_nwis <- eligible_sites %>%
  str_subset("^USGS-") %>%
  str_remove("^USGS-")

batch_size <- 50
batches    <- split(usgs_sites_nwis, ceiling(seq_along(usgs_sites_nwis) / batch_size))

nwis_list <- map(seq_along(batches), function(i) {
  batch <- batches[[i]]
  cat("  NWIS batch", i, "of", length(batches), "(", length(batch), "sites)...\n")

  result <- tryCatch(
    readNWISdv(
      siteNumbers = batch,
      parameterCd = c("00060", "00010"),
      startDate   = START_DATE,
      endDate     = END_DATE
    ),
    error = function(e) { cat("  NWIS ERROR:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(result) || nrow(result) == 0) {
    cat("  Batch", i, "no data\n")
    return(NULL)
  }

  renameNWISColumns(result) %>%
    select(site_no, Date, any_of(c("Flow", "Wtemp"))) %>%
    mutate(
      Date    = as.Date(Date),
      Flow    = suppressWarnings(as.numeric(Flow)),
      Wtemp   = suppressWarnings(as.numeric(Wtemp)),
      site_no = str_pad(as.character(site_no), 8, pad = "0")
    ) %>%
    filter(!is.na(Flow) | !is.na(Wtemp))
})

nwis_clean <- bind_rows(nwis_list)
cat("Total NWIS rows:", nrow(nwis_clean), "\n")
saveRDS(nwis_clean, file.path(OUT_DIR, "nwis_daily_slim.rds"))
