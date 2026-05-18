# =============================================================================
# Aggregate SINAN Cases + Population -> Municipality-Year Rates
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Outputs:
#   agg_all     : all municipality-years (incl. zeros) — for mapping & LISA
#   agg_nonzero : municipality-years with at least 1 case — for descriptives
#
# Key decisions:
#   - Unit of analysis: municipality of RESIDENCE (muni_resid / v19)
#     Rationale: captures where transmission occurred, not where reported
#   - Denominator: annual population from population file (POPULACAO)
#   - Incidence rate: cases per 100,000 inhabitants
#   - CFR: VL deaths / confirmed cases  (2007-2018 only; NA for 2001-2006)
#   - Excluded: 33 unmatched codes (state-level, missing, unknown)
#
# =============================================================================

library(dplyr)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"

# =============================================================================
# 1. LOAD
# =============================================================================

sinan <- readRDS(file.path(proc_dir, "sinan_confirmed_2001_2018.rds"))
pop   <- readRDS(file.path(proc_dir, "population_2001_2018.rds"))

message("SINAN confirmed: ", nrow(sinan), " cases")
message("Population rows: ", nrow(pop))

# =============================================================================
# 2. IDENTIFY VALID MUNICIPALITY CODES
# =============================================================================
# Keep only codes that exist in the population file — this excludes all 33
# unmatched codes (state-level, missing, unknown).

valid_codes <- unique(pop$MUNIC_RES)

sinan_clean <- sinan |>
  filter(muni_resid %in% valid_codes)

n_excluded <- nrow(sinan) - nrow(sinan_clean)
message("\nCases excluded (unmatched municipality code): ", n_excluded,
        " (", round(100 * n_excluded / nrow(sinan), 2), "%)")
message("Cases retained: ", nrow(sinan_clean))

# =============================================================================
# 3. AGGREGATE CASES BY MUNICIPALITY-YEAR
# =============================================================================

cases_agg <- sinan_clean |>
  group_by(muni_resid, year) |>
  summarise(
    n_cases    = n(),
    n_vl_death = sum(vl_death, na.rm = TRUE),   # NA for 2001-2006 rows
    n_hiv      = sum(hiv == 1, na.rm = TRUE),
    n_male     = sum(sex == "M", na.rm = TRUE),
    n_female   = sum(sex == "F", na.rm = TRUE),
    n_under5   = sum(age_years < 5, na.rm = TRUE),
    n_over40   = sum(age_years >= 40, na.rm = TRUE),
    has_death_data = any(period == "2007-2018"),  # TRUE if CFR computable
    .groups = "drop"
  )

message("\nMunicipality-year combinations with >= 1 case: ", nrow(cases_agg))

# =============================================================================
# 4. BUILD COMPLETE MUNICIPALITY-YEAR GRID (include zeros)
# =============================================================================
# For spatial analysis (LISA, CAR models) every municipality needs a value
# in every year, even if cases = 0.

all_munis <- unique(pop$MUNIC_RES)
all_years <- 2001:2018

grid <- expand.grid(
  muni_resid = all_munis,
  year       = all_years,
  stringsAsFactors = FALSE
)

message("Full municipality-year grid: ", nrow(grid), " rows")

# =============================================================================
# 5. JOIN CASES + POPULATION ONTO GRID
# =============================================================================

agg <- grid |>
  left_join(cases_agg,
            by = c("muni_resid", "year")) |>
  left_join(pop,
            by = c("muni_resid" = "MUNIC_RES", "year" = "ANO")) |>
  mutate(
    # Fill zero cases for municipalities with no reports that year
    n_cases    = coalesce(n_cases,    0L),
    n_vl_death = coalesce(n_vl_death, 0L),
    n_hiv      = coalesce(n_hiv,      0L),
    n_male     = coalesce(n_male,     0L),
    n_female   = coalesce(n_female,   0L),
    n_under5   = coalesce(n_under5,   0L),
    n_over40   = coalesce(n_over40,   0L),
    
    # Incidence rate per 100,000
    incidence_rate = ifelse(
      !is.na(POPULACAO) & POPULACAO > 0,
      n_cases / POPULACAO * 100000,
      NA_real_
    ),
    
    # Case fatality rate (proportion) — only for 2007-2018 and where cases > 0
    cfr = case_when(
      year < 2007        ~ NA_real_,   # death code not available
      n_cases == 0       ~ NA_real_,   # undefined
      TRUE               ~ n_vl_death / n_cases
    ),
    
    # State code (first 2 digits of municipality code)
    state_code = substr(muni_resid, 1, 2),
    
    # Period flag for period-level analyses
    period = case_when(
      year <= 2006 ~ "2001-2006",
      year <= 2012 ~ "2007-2012",
      TRUE         ~ "2013-2018"
    )
  )

message("\nFull aggregated dataset: ", nrow(agg), " municipality-year rows")
message("Municipalities with missing population: ",
        sum(is.na(agg$POPULACAO)))

# =============================================================================
# 6. QUALITY CHECKS
# =============================================================================

message("\n--- Cases per year (should match SINAN) ---")
print(
  agg |>
    group_by(year) |>
    summarise(total_cases = sum(n_cases), .groups = "drop"),
  n = 18
)

message("\n--- National incidence rate per year (cases/100k) ---")
nat_rates <- agg |>
  group_by(year) |>
  summarise(
    total_cases = sum(n_cases),
    total_pop   = sum(POPULACAO, na.rm = TRUE),
    nat_rate    = total_cases / total_pop * 100000,
    .groups = "drop"
  )
print(nat_rates, n = 18)

message("\n--- VL deaths per year (2007-2018) ---")
print(
  agg |>
    filter(year >= 2007) |>
    group_by(year) |>
    summarise(total_deaths = sum(n_vl_death), .groups = "drop"),
  n = 12
)

message("\n--- Incidence rate distribution (non-zero municipalities) ---")
nonzero <- agg |> filter(n_cases > 0)
message("  Municipality-years with >= 1 case: ", nrow(nonzero))
message("  Mean incidence rate:   ", round(mean(nonzero$incidence_rate, na.rm=TRUE), 2))
message("  Median incidence rate: ", round(median(nonzero$incidence_rate, na.rm=TRUE), 2))
message("  Max incidence rate:    ", round(max(nonzero$incidence_rate, na.rm=TRUE), 2))

# =============================================================================
# 7. SAVE
# =============================================================================

# Full grid (zeros included) — for LISA, CAR models, maps
saveRDS(agg, file.path(proc_dir, "agg_municipality_year.rds"))

# Non-zero only — for descriptive tables and CFR analyses
agg_nonzero <- agg |> filter(n_cases > 0)
saveRDS(agg_nonzero, file.path(proc_dir, "agg_nonzero.rds"))

message("\nSaved:")
message("  agg_municipality_year.rds  (", nrow(agg), " rows — includes zeros)")
message("  agg_nonzero.rds            (", nrow(agg_nonzero), " rows — cases > 0)")
message("\nNext: run 04_bayesian_smoothing.R")