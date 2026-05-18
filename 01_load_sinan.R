# =============================================================================
# Load and Combine SINAN VL Case Data — 2001–2006 and 2007–2018
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Data locations (confirmed from screenshots):
#   D:\AIMS\Research phase\R\data\DISTRIB_VI_LEISHMANIASIS\
#     VL_Brazil_2001_2006_Ali_edited.csv
#     VL_Brazil_2007_2018_Ali_edited.csv
#
# Variable encoding confirmed by inspection:
#
#   v6  : year of notification
#   v9  : municipality of notification (5 digits)
#   v19 : municipality of residence    (5 digits)
#   v12 : age
#         2001-2006 — letter prefix: A=years  M=months  D=days  I=ignored
#                     e.g. A004=4yrs  M006=6mo  D023=23days
#         2007-2018 — digit prefix:  4=years  3=months  2=days  1=hours
#                     e.g. 4006=6yrs  3003=3mo  2009=9days
#   v13 : sex          M / F / I  (identical both files)
#   v32 : HIV          1=yes  2=no  9=ignored  (identical both files)
#   v43 : final classification  1=confirmed  2=discarded  (identical both files)
#   v48 : outcome
#         2001-2006 — 1=cure  2=default  9=ignored  (NO VL-death code)
#         2007-2018 — 1=cure  2=default  3=VL death  4=other-cause death  5=transferred
#
#   Municipality code join key:
#     SINAN v9/v19 = 5 digits (e.g. 24081)
#     Population MUNIC_RES = 6 digits (e.g. 240810)
#     Shapefile code6 = 6 digits (first 6 of 7-digit IBGE code)
#     Conversion: multiply SINAN 5-digit by 10 to get 6-digit
#                 e.g. 24081 × 10 = 240810  ✓
#
# =============================================================================

library(dplyr)
library(readr)

# =============================================================================
# 1. PATHS — adjust only if your R working directory is not the project root
# =============================================================================

data_dir <- "D:/AIMS/Research phase/R/data/DISTRIB_VI_LEISHMANIASIS"

path_0106 <- file.path(data_dir, "VL_Brazil_2001_2006_Ali_edited.csv")
path_0718 <- file.path(data_dir, "VL_Brazil_2007_2018_Ali_edited.csv")

# Verify files exist before proceeding
stopifnot(
  "2001-2006 CSV not found — check data_dir path" = file.exists(path_0106),
  "2007-2018 CSV not found — check data_dir path" = file.exists(path_0718)
)
message("Both CSV files found. Reading...")

# =============================================================================
# 2. AGE PARSING FUNCTIONS
# =============================================================================

# 2001-2006: letter prefix (A/M/D/I) + 3-digit number
parse_age_0106 <- function(v12) {
  unit   <- substr(trimws(v12), 1, 1)
  number <- suppressWarnings(as.integer(substr(trimws(v12), 2, 4)))
  case_when(
    unit == "A" ~ number,
    unit == "M" ~ floor(number / 12),
    unit == "D" ~ 0L,
    TRUE         ~ NA_integer_
  )
}

# 2007-2018: digit prefix (4/3/2/1) + 3-digit number
parse_age_0718 <- function(v12) {
  unit   <- substr(trimws(v12), 1, 1)
  number <- suppressWarnings(as.integer(substr(trimws(v12), 2, 4)))
  case_when(
    unit == "4" ~ number,
    unit == "3" ~ floor(number / 12),
    unit == "2" ~ 0L,
    unit == "1" ~ 0L,
    TRUE         ~ NA_integer_
  )
}

# 5-digit SINAN code -> 6-digit code matching population file and shapefile
to_code6 <- function(x) {
  formatC(as.integer(trimws(x)) * 10L, width = 6, flag = "0")
}

# =============================================================================
# 3. READ 2001–2006
# =============================================================================

message("Reading 2001-2006...")

raw_0106 <- read_csv(
  path_0106,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

message("  Raw rows: ", nrow(raw_0106))  # expect 33,479

sinan_0106 <- raw_0106 |>
  transmute(
    year        = as.integer(trimws(v6)),
    muni_notif  = to_code6(v9),
    muni_resid  = to_code6(v19),
    age_years   = parse_age_0106(v12),
    sex         = case_when(
      trimws(v13) == "M" ~ "M",
      trimws(v13) == "F" ~ "F",
      TRUE               ~ NA_character_
    ),
    hiv         = case_when(
      trimws(v32) == "1" ~ 1L,
      trimws(v32) == "2" ~ 0L,
      TRUE               ~ NA_integer_
    ),
    confirmed   = as.integer(trimws(v43) == "1"),
    vl_death    = NA_integer_,        # outcome code unavailable 2001-2006
    outcome_raw = trimws(v48),
    period      = "2001-2006"
  )

message("  Parsed rows: ", nrow(sinan_0106))

# =============================================================================
# 4. READ 2007–2018
# =============================================================================

message("Reading 2007-2018...")

raw_0718 <- read_csv(
  path_0718,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

message("  Raw rows: ", nrow(raw_0718))  # expect 112,962

sinan_0718 <- raw_0718 |>
  transmute(
    year        = as.integer(trimws(v6)),
    muni_notif  = to_code6(v9),
    muni_resid  = to_code6(v19),
    age_years   = parse_age_0718(v12),
    sex         = case_when(
      trimws(v13) == "M" ~ "M",
      trimws(v13) == "F" ~ "F",
      TRUE               ~ NA_character_
    ),
    hiv         = case_when(
      trimws(v32) == "1" ~ 1L,
      trimws(v32) == "2" ~ 0L,
      TRUE               ~ NA_integer_
    ),
    confirmed   = as.integer(trimws(v43) == "1"),
    vl_death    = as.integer(trimws(v48) == "3"),
    outcome_raw = trimws(v48),
    period      = "2007-2018"
  )

message("  Parsed rows: ", nrow(sinan_0718))

# =============================================================================
# 5. COMBINE
# =============================================================================

sinan_all <- bind_rows(sinan_0106, sinan_0718)

message("\n=== Combined dataset ===")
message("Total rows:       ", nrow(sinan_all))
message("Year range:       ", min(sinan_all$year, na.rm = TRUE),
        " - ", max(sinan_all$year, na.rm = TRUE))
message("Confirmed cases:  ", sum(sinan_all$confirmed, na.rm = TRUE))
message("VL deaths (2007-2018 only): ",
        sum(sinan_all$vl_death, na.rm = TRUE))

# =============================================================================
# 6. AGE GROUPS
# =============================================================================

sinan_all <- sinan_all |>
  mutate(
    age_group = case_when(
      age_years <  1                    ~ "<1",
      age_years >= 1  & age_years <  5  ~ "1-4",
      age_years >= 5  & age_years < 10  ~ "5-9",
      age_years >= 10 & age_years < 20  ~ "10-19",
      age_years >= 20 & age_years < 40  ~ "20-39",
      age_years >= 40 & age_years < 60  ~ "40-59",
      age_years >= 60                   ~ "60+",
      TRUE                              ~ NA_character_
    )
  )

# =============================================================================
# 7. CONFIRMED-ONLY SUBSET
# =============================================================================

sinan_confirmed <- sinan_all |>
  filter(confirmed == 1)

message("\nConfirmed cases only: ", nrow(sinan_confirmed))

# =============================================================================
# 8. QUALITY CHECKS
# =============================================================================

message("\n--- Year distribution (confirmed) ---")
print(table(sinan_confirmed$year))

message("\n--- Sex (confirmed) ---")
print(table(sinan_confirmed$sex, useNA = "ifany"))

message("\n--- HIV status (confirmed) ---")
print(table(sinan_confirmed$hiv, useNA = "ifany"))

message("\n--- Age group (confirmed) ---")
print(table(sinan_confirmed$age_group, useNA = "ifany"))

message("\n--- VL deaths (confirmed, 2007-2018 only) ---")
print(table(
  sinan_confirmed$vl_death[sinan_confirmed$period == "2007-2018"],
  useNA = "ifany"
))

# =============================================================================
# 9. COVERAGE CHECK AGAINST POPULATION DATA
# =============================================================================
# Uncomment after loading population data (next script).
# Verifies that the 5->6 digit conversion is correct.
#
# pop_codes   <- unique(pop$MUNIC_RES)           # already 6-digit
# sinan_codes <- unique(sinan_confirmed$muni_resid)
# missing     <- setdiff(sinan_codes, pop_codes)
# message("SINAN codes not matched in population file: ", length(missing))
# if (length(missing) > 0) print(head(missing, 20))

# =============================================================================
# 10. SAVE
# =============================================================================

out_dir <- "D:/AIMS/Research phase/R/data/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(sinan_all,
        file.path(out_dir, "sinan_all_2001_2018.rds"))
saveRDS(sinan_confirmed,
        file.path(out_dir, "sinan_confirmed_2001_2018.rds"))

message("\nSaved to: ", out_dir)
message("  sinan_all_2001_2018.rds       (", nrow(sinan_all), " rows)")
message("  sinan_confirmed_2001_2018.rds (", nrow(sinan_confirmed), " rows)")
message("\nNext: run 02_load_population.R")