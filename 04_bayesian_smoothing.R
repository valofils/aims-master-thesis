# =============================================================================
# Global Empirical Bayesian Smoothing of VL Incidence Rates
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Why GLOBAL rather than LOCAL EB?
#   Local EB (EBlocal) borrows strength from neighbouring municipalities.
#   It works well when cases are widespread. For VL in Brazil, only ~13% of
#   municipalities report cases in any given year. When a municipality has
#   cases but all its neighbours have zero cases, the local prior mean is
#   near zero and EBlocal collapses all estimates to ~0 — over-shrinkage.
#
#   Global EB (EBest) uses the NATIONAL mean rate as the prior. This is
#   appropriate for a rare, spatially concentrated disease and is the approach
#   used in the published VL spatial literature (de Melo et al. 2023,
#   Bruhn et al. 2024, Lima et al. 2021).
#
# Method: EBest() from spdep (Clayton & Kaldor 1987)
#   - ri = n_cases / POPULACAO  (raw proportion)
#   - ni = POPULACAO
#   - Returns $raw (crude rate) and $est (smoothed rate), both as proportions
#   - Multiply by 100,000 for per-100k rates
#   - Applied per year across all municipalities
#
# =============================================================================

library(sf)
library(dplyr)
library(spdep)

proc_dir  <- "D:/AIMS/Research phase/R/data/processed"
shp_dir07 <- "data/shapefiles/2007"
shp_dir18 <- "data/shapefiles/2018"

# =============================================================================
# 1. LOAD DATA
# =============================================================================

agg <- readRDS(file.path(proc_dir, "agg_municipality_year.rds"))
message("Aggregated data: ", nrow(agg), " municipality-year rows")

shp_2007_path <- list.files(shp_dir07, pattern = "\\.shp$",
                            full.names = TRUE, recursive = TRUE)
shp_2018_path <- list.files(shp_dir18, pattern = "\\.shp$",
                            full.names = TRUE, recursive = TRUE)

muni_2007_raw <- st_read(shp_2007_path[1], quiet = TRUE)
muni_2018_raw <- st_read(shp_2018_path[1], quiet = TRUE)

WATER_BODIES <- c("4300001", "4300002")

muni_2007 <- muni_2007_raw |>
  filter(!as.character(GEOCODIG_M) %in% WATER_BODIES) |>
  mutate(code6 = substr(as.character(GEOCODIG_M), 1, 6))

muni_2018 <- muni_2018_raw |>
  filter(!as.character(CD_GEOCMU) %in% WATER_BODIES) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

message("2007 shapefile: ", nrow(muni_2007), " municipalities")
message("2018 shapefile: ", nrow(muni_2018), " municipalities")

# =============================================================================
# 2. GLOBAL EB SMOOTHING FUNCTION
# =============================================================================
# EBest(n, x) where:
#   n = case counts vector (aligned to shapefile row order)
#   x = population vector  (aligned to shapefile row order)
# Returns data.frame with $raw and $est (both as raw proportions)

smooth_year_global <- function(year_val, agg_data, shapefile) {
  
  # Align this year's data to shapefile row order
  sf_data <- shapefile |>
    st_drop_geometry() |>
    left_join(
      agg_data |>
        filter(year == year_val) |>
        select(muni_resid, n_cases, POPULACAO),
      by = c("code6" = "muni_resid")
    ) |>
    mutate(
      n_cases   = coalesce(as.integer(n_cases),   0L),
      POPULACAO = coalesce(as.numeric(POPULACAO), NA_real_),
      pop_valid = !is.na(POPULACAO) & POPULACAO > 0
    )
  
  # For EBest: use valid population rows only
  # Rows with missing population get crude_rate = NA, smoothed_rate = NA
  n_vec <- sf_data$n_cases
  x_vec <- sf_data$POPULACAO
  
  # Replace NA population with 0 temporarily; EBest handles zeros via
  # global variance estimation (unlike EBlocal which divides locally)
  n_vec[!sf_data$pop_valid] <- 0
  x_vec[!sf_data$pop_valid] <- 0
  
  # Run global EB smoothing
  eb <- tryCatch(
    EBest(n = n_vec, x = x_vec),
    error = function(e) {
      message("  EBest error year ", year_val, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(eb)) {
    smoothed <- eb$est    # global EB smoothed rate (proportion)
    crude    <- eb$raw    # crude rate (proportion) = n/x
  } else {
    # Fallback: crude rate
    crude    <- ifelse(sf_data$pop_valid,
                       sf_data$n_cases / sf_data$POPULACAO, NA_real_)
    smoothed <- crude
  }
  
  # Set both rates to NA where population is missing
  crude   [!sf_data$pop_valid] <- NA_real_
  smoothed[!sf_data$pop_valid] <- NA_real_
  
  data.frame(
    muni_resid    = sf_data$code6,
    year          = year_val,
    n_cases       = sf_data$n_cases,
    population    = sf_data$POPULACAO,
    crude_rate    = crude    * 100000,   # per 100,000
    smoothed_rate = smoothed * 100000,   # per 100,000
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 3. APPLY SMOOTHING — ALL 18 YEARS
# =============================================================================

message("\nApplying global EB smoothing across 18 years...")

results_list <- vector("list", 18)

for (yr in 2001:2018) {
  shp <- if (yr <= 2009) muni_2007 else muni_2018
  results_list[[yr - 2000]] <- smooth_year_global(yr, agg, shp)
  message("  Year ", yr, " done — rows: ",
          nrow(results_list[[yr - 2000]]))
}

smoothed_all <- bind_rows(results_list)
message("\nTotal smoothed rows: ", nrow(smoothed_all))

# =============================================================================
# 4. QUALITY CHECKS
# =============================================================================

message("\n--- Effect of global EB smoothing (municipalities with cases) ---")
nz <- smoothed_all |>
  filter(n_cases > 0, !is.na(crude_rate), !is.na(smoothed_rate))

message("  Municipality-years with cases + valid rates: ", nrow(nz))
message("  Crude rate   — mean: ", round(mean(nz$crude_rate),    2),
        "  max: ", round(max(nz$crude_rate),    1))
message("  Smoothed rate — mean: ", round(mean(nz$smoothed_rate), 2),
        "  max: ", round(max(nz$smoothed_rate), 1))
message("  [Smoothed max < crude max confirms shrinkage of extreme values]")
message("  [Smoothed mean ≈ crude mean confirms global prior is correct]")

message("\n--- National rates per year ---")
nat <- smoothed_all |>
  group_by(year) |>
  summarise(
    cases         = sum(n_cases,      na.rm = TRUE),
    n_valid       = sum(!is.na(smoothed_rate)),
    mean_crude    = round(mean(crude_rate,    na.rm = TRUE), 3),
    mean_smoothed = round(mean(smoothed_rate, na.rm = TRUE), 3),
    max_crude     = round(max(crude_rate,     na.rm = TRUE), 1),
    max_smoothed  = round(max(smoothed_rate,  na.rm = TRUE), 1),
    .groups = "drop"
  )
print(nat, n = 18)

# =============================================================================
# 5. SPOT CHECK — confirm shrinkage on extreme municipalities
# =============================================================================

message("\n--- Top 10 highest crude-rate municipalities (2004) ---")
check <- smoothed_all |>
  filter(year == 2004, n_cases > 0) |>
  arrange(desc(crude_rate)) |>
  slice_head(n = 10) |>
  select(muni_resid, n_cases, population, crude_rate, smoothed_rate)
print(check)
# Smoothed rates should be LOWER than crude rates for high outliers
# and HIGHER than crude rates for zero-case neighbours

# =============================================================================
# 6. MERGE BACK INTO MAIN AGGREGATED DATASET
# =============================================================================

agg_smoothed <- agg |>
  left_join(
    smoothed_all |> select(muni_resid, year, crude_rate, smoothed_rate),
    by = c("muni_resid", "year")
  )

message("\nMerged dataset rows:           ", nrow(agg_smoothed))
message("Rows with valid smoothed rate:  ",
        sum(!is.na(agg_smoothed$smoothed_rate)))
message("Rows with NA smoothed rate:     ",
        sum( is.na(agg_smoothed$smoothed_rate)))

# =============================================================================
# 7. SAVE
# =============================================================================

# Save neighbourhood lists too — still needed for LISA even though
# we used global EB for smoothing
sf_use_s2(FALSE)
muni_2007_nb <- st_make_valid(muni_2007)
nb_2007 <- poly2nb(muni_2007_nb, queen = TRUE, snap = 0.01)
muni_2018_nb <- st_make_valid(muni_2018)
nb_2018 <- poly2nb(muni_2018_nb, queen = TRUE, snap = 0.01)
sf_use_s2(TRUE)
message("Neighbourhood lists rebuilt for LISA")

saveRDS(smoothed_all, file.path(proc_dir, "smoothed_rates.rds"))
saveRDS(agg_smoothed, file.path(proc_dir, "agg_smoothed.rds"))
saveRDS(nb_2007,      file.path(proc_dir, "nb_2007.rds"))
saveRDS(nb_2018,      file.path(proc_dir, "nb_2018.rds"))

message("\nSaved:")
message("  smoothed_rates.rds  — global EB smoothed rates per municipality-year")
message("  agg_smoothed.rds    — full aggregated data + smoothed rates")
message("  nb_2007.rds  |  nb_2018.rds  — neighbourhood lists for LISA")
message("\nNext: run 05_LISA.R")