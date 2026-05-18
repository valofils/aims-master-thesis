# =============================================================================
# Load Population Data — Brazil Municipalities 2001–2018
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# File formats confirmed by inspection:
#   POPTBR01–13 : ZIP containing CSV  (MUNIC_RES = 6 digits)
#   POPTBR14–16 : ZIP containing DBF  (MUNIC_RES = 6 digits)
#   POPTBR17–18 : ZIP containing DBF  (MUNIC_RES = 7 digits — truncate to 6)
#
# Requires package: foreign (reads DBF — ships with base R, no install needed)
#
# =============================================================================

library(dplyr)
library(readr)
library(foreign)   # read.dbf() — included in base R

pop_dir <- "D:/AIMS/Research phase/R/data/Pop_Brazil_2001_2018"

# =============================================================================
# HELPER: normalise MUNIC_RES to 6 digits
# =============================================================================

norm_code6 <- function(x) {
  x <- trimws(as.character(x))
  nc <- nchar(x)
  case_when(
    nc == 6 ~ x,                                    # already 6 digits
    nc == 7 ~ substr(x, 1, 6),                      # 7-digit: take first 6
    nc  < 6 ~ formatC(as.integer(x), width = 6,    # short: left-pad with 0
                      flag = "0"),
    TRUE    ~ NA_character_
  )
}

# =============================================================================
# LOAD ALL 18 YEARS
# =============================================================================

pop_list <- vector("list", 18)

for (i in seq_len(18)) {
  
  year     <- 2000L + i
  zip_name <- sprintf("POPTBR%02d.zip", i)
  zip_path <- file.path(pop_dir, zip_name)
  
  if (!file.exists(zip_path)) {
    warning("Not found, skipping: ", zip_path)
    next
  }
  
  # Peek inside the ZIP to decide format
  zip_contents <- unzip(zip_path, list = TRUE)$Name
  has_csv <- any(grepl("\\.csv$", zip_contents, ignore.case = TRUE))
  has_dbf <- any(grepl("\\.dbf$", zip_contents, ignore.case = TRUE))
  
  tmp <- tempfile()
  dir.create(tmp)
  
  if (has_csv) {
    # ── CSV format (years 2001–2013) ─────────────────────────────────────────
    csv_name <- zip_contents[grepl("\\.csv$", zip_contents, ignore.case = TRUE)][1]
    df <- read_csv(
      unz(zip_path, csv_name),
      col_types = cols(.default = "c"),
      show_col_types = FALSE
    )
    names(df) <- toupper(names(df))
    df <- df |>
      transmute(
        MUNIC_RES = norm_code6(MUNIC_RES),
        ANO       = as.integer(ANO),
        POPULACAO = as.numeric(POPULACAO)
      )
    
  } else if (has_dbf) {
    # ── DBF format (years 2014–2018) ─────────────────────────────────────────
    dbf_name <- zip_contents[grepl("\\.dbf$", zip_contents, ignore.case = TRUE)][1]
    unzip(zip_path, files = dbf_name, exdir = tmp, overwrite = TRUE)
    dbf_path <- file.path(tmp, dbf_name)
    
    raw <- read.dbf(dbf_path, as.is = TRUE)
    names(raw) <- toupper(names(raw))
    df <- raw |>
      as_tibble() |>
      transmute(
        MUNIC_RES = norm_code6(MUNIC_RES),
        ANO       = as.integer(ANO),
        POPULACAO = as.numeric(POPULACAO)
      )
    
  } else {
    warning("No CSV or DBF found in ", zip_name)
    next
  }
  
  # Remove NA codes and zero populations
  df <- df |>
    filter(!is.na(MUNIC_RES), !is.na(POPULACAO), POPULACAO > 0)
  
  pop_list[[i]] <- df
  message("Loaded: ", zip_name,
          "  year=", year,
          "  rows=", nrow(df),
          "  format=", ifelse(has_csv, "CSV", "DBF"))
  
  unlink(tmp, recursive = TRUE)
}

# =============================================================================
# STACK
# =============================================================================

pop <- bind_rows(pop_list)

message("\n=== Population data combined ===")
message("Total rows:    ", nrow(pop))
message("Year range:    ", min(pop$ANO), " - ", max(pop$ANO))
message("Unique munis:  ", n_distinct(pop$MUNIC_RES))

# =============================================================================
# QUALITY CHECKS
# =============================================================================

message("\n--- Rows per year ---")
print(table(pop$ANO))

message("\n--- Population range ---")
message("  Min:  ", format(min(pop$POPULACAO),  big.mark = ","))
message("  Max:  ", format(max(pop$POPULACAO),  big.mark = ","))
message("  Mean: ", format(round(mean(pop$POPULACAO)), big.mark = ","))

message("\n--- Year coverage per municipality ---")
year_counts <- pop |>
  group_by(MUNIC_RES) |>
  summarise(n_years = n_distinct(ANO), .groups = "drop")
print(table(year_counts$n_years))

# =============================================================================
# COVERAGE CHECK AGAINST SINAN
# =============================================================================

sinan_confirmed <- readRDS(
  "D:/AIMS/Research phase/R/data/processed/sinan_confirmed_2001_2018.rds"
)

sinan_codes      <- unique(sinan_confirmed$muni_resid)
pop_codes        <- unique(pop$MUNIC_RES)
missing_in_pop   <- setdiff(sinan_codes, pop_codes)
missing_in_sinan <- setdiff(pop_codes, sinan_codes)

message("\n--- Code coverage ---")
message("SINAN unique muni_resid codes:     ", length(sinan_codes))
message("Population unique MUNIC_RES codes: ", length(pop_codes))
message("SINAN codes NOT in population:     ", length(missing_in_pop))
message("Population codes NOT in SINAN:     ", length(missing_in_sinan))

if (length(missing_in_pop) > 0) {
  message("  Unmatched SINAN codes:")
  print(sort(missing_in_pop))
}

# =============================================================================
# SAVE
# =============================================================================

out_dir <- "D:/AIMS/Research phase/R/data/processed"
saveRDS(pop, file.path(out_dir, "population_2001_2018.rds"))

message("\nSaved: ", file.path(out_dir, "population_2001_2018.rds"))
message("  Rows: ", nrow(pop))
message("\nNext: run 03_aggregate.R")