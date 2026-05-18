# =============================================================================
# Local Indicators of Spatial Association (LISA) — VL Case Fatality Rate
# CORRECTED: Pooled 2007–2018 (matching INLA RQ2 structure)
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# CORRECTION RATIONALE:
# ====================
# Your INLA RQ2 models (Model A & B) pool CFR data across 2007–2018,
# computing a SINGLE municipality-level CFR estimate for the entire period.
# 
# For consistency, LISA should also use pooled 2007–2018 CFR data
# (not split into 2007–2012 and 2013–2018 periods).
#
# This ensures:
# 1. LISA and INLA use identical CFR denominators
# 2. Hotspots/coldspots are identified in the same data as INLA spatial effects
# 3. Direct comparison between LISA and INLA results
#
# Cluster types:
#   High-High (HH): hotspot — high CFR surrounded by high CFR neighbours
#   Low-Low   (LL): coldspot — low CFR surrounded by low CFR neighbours
#   High-Low  (HL): spatial outlier — isolated high CFR
#   Low-High  (LH): spatial outlier — isolated low CFR
#   Not significant: p >= 0.05
#
# spdep version: 1.4.2
#   localmoran_perm() with nsim=999 for permutation-based p-values
#
# =============================================================================

library(sf)
library(dplyr)
library(spdep)
library(tidyr)

proc_dir  <- "D:/AIMS/Research phase/R/data/processed"

set.seed(2024)

# =============================================================================
# 1. LOAD DATA AND SHAPEFILES
# =============================================================================

agg <- readRDS(file.path(proc_dir, "agg_municipality_year.rds"))
nb_2018 <- readRDS(file.path(proc_dir, "nb_2018.rds"))

message("Aggregated data: ", nrow(agg), " municipality-year rows")

shp_dir18 <- "data/shapefiles/2018"
shp_2018_path <- list.files(shp_dir18, pattern = "\\.shp$",
                             full.names = TRUE, recursive = TRUE)

muni_2018_raw <- st_read(shp_2018_path[1], quiet = TRUE)

# Exclude water bodies
WATER_BODIES <- c("4300001", "4300002")

muni_2018 <- muni_2018_raw |>
  filter(!as.character(CD_GEOCMU) %in% WATER_BODIES) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

message("2018 shapefile: ", nrow(muni_2018), " municipalities")

# =============================================================================
# 2. AGGREGATE CFR ACROSS ENTIRE 2007–2018 PERIOD
# =============================================================================
#
# Following your INLA RQ2 approach: pool CFR across all years 2007–2018
# to get one CFR estimate per municipality.
#
# This differs from RQ1, which has annual data (2001–2018).
# RQ2 CFR has no temporal variation within the model — it's a static
# municipality-level outcome.

message("\n=== CFR Data Summary (2007–2018 pooled) ===")

cfr_data_pooled <- agg |>
  filter(year >= 2007) |>
  group_by(muni_resid) |>
  summarise(
    total_cases = sum(n_cases, na.rm = TRUE),
    total_deaths = sum(n_vl_death, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cfr = ifelse(total_cases > 0, total_deaths / total_cases, NA_real_),
    cfr_pct = cfr * 100
  )

message("Valid CFR observations (municipalities with >= 1 case): ",
        sum(!is.na(cfr_data_pooled$cfr)))
message("CFR range: ",
        round(min(cfr_data_pooled$cfr, na.rm = TRUE), 4), " to ",
        round(max(cfr_data_pooled$cfr, na.rm = TRUE), 4))
message("CFR pct range: ",
        round(min(cfr_data_pooled$cfr_pct, na.rm = TRUE), 2), "% to ",
        round(max(cfr_data_pooled$cfr_pct, na.rm = TRUE), 2), "%")

message("\nCases distribution:")
print(
  cfr_data_pooled |>
    summarise(
      n_municipalities = n(),
      n_with_cases = sum(total_cases > 0),
      total_cases_national = sum(total_cases),
      total_deaths_national = sum(total_deaths),
      national_cfr = total_deaths_national / total_cases_national
    )
)

# =============================================================================
# 3. CORE LISA FUNCTION — spdep 1.4.2 compatible
# =============================================================================

run_lisa_cfr <- function(cfr_vec, nb, shapefile, label = "") {
  
  n <- length(cfr_vec)
  stopifnot(length(nb) == n)
  
  # Replace NA with 0 — municipalities with zero cases cannot form CFR clusters
  cfr_clean <- ifelse(is.na(cfr_vec), 0, cfr_vec)
  
  # Row-standardised spatial weights
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # localmoran_perm: permutation-based p-values (spdep >= 1.2)
  local_moran <- tryCatch(
    localmoran_perm(
      x           = cfr_clean,
      listw       = lw,
      zero.policy = TRUE,
      alternative = "two.sided",
      nsim        = 999,
      iseed       = 2024
    ),
    error = function(e) {
      message("  Falling back to localmoran() (analytical): ", e$message)
      localmoran(
        x           = cfr_clean,
        listw       = lw,
        zero.policy = TRUE,
        alternative = "two.sided"
      )
    }
  )
  
  local_i <- local_moran[, "Ii"]
  
  # p-value column name differs between localmoran and localmoran_perm
  p_col   <- grep("^Pr", colnames(local_moran), value = TRUE)[1]
  p_value <- local_moran[, p_col]
  
  # Classify quadrant
  x_mean <- mean(cfr_clean, na.rm = TRUE)
  lag_x  <- lag.listw(lw, cfr_clean, zero.policy = TRUE)
  
  cluster_type <- case_when(
    p_value >= 0.05                        ~ "Not significant",
    cfr_clean >= x_mean & lag_x >= x_mean ~ "High-High",
    cfr_clean <  x_mean & lag_x <  x_mean ~ "Low-Low",
    cfr_clean >= x_mean & lag_x <  x_mean ~ "High-Low",
    cfr_clean <  x_mean & lag_x >= x_mean ~ "Low-High",
    TRUE                                   ~ "Not significant"
  )
  
  if (label != "")
    message("  ", label,
            "  HH=", sum(cluster_type == "High-High"),
            "  LL=", sum(cluster_type == "Low-Low"),
            "  HL=", sum(cluster_type == "High-Low"),
            "  LH=", sum(cluster_type == "Low-High"),
            "  NS=", sum(cluster_type == "Not significant"))
  
  data.frame(
    muni_resid   = shapefile$code6,
    local_i      = local_i,
    p_value      = p_value,
    lag_cfr      = lag_x,
    cluster_type = cluster_type,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 4. RUN LISA ON POOLED 2007–2018 CFR
# =============================================================================

message("\n=== LISA Analysis (CFR, 2007–2018 pooled) ===\n")

# Align CFR data to shapefile row order
shp_data <- muni_2018 |>
  st_drop_geometry() |>
  left_join(cfr_data_pooled, by = c("code6" = "muni_resid"))

# Run LISA
lisa_cfr <- run_lisa_cfr(shp_data$cfr, nb_2018, muni_2018,
                         label = "2007-2018 pooled")
lisa_cfr$total_cases <- shp_data$total_cases
lisa_cfr$total_deaths <- shp_data$total_deaths
lisa_cfr$cfr <- shp_data$cfr
lisa_cfr$cfr_pct <- shp_data$cfr_pct

# Summary
message("\n--- LISA cluster counts (2007–2018 pooled) ---")
print(
  lisa_cfr |>
    group_by(cluster_type) |>
    summarise(n = n(), .groups = "drop") |>
    arrange(factor(cluster_type, 
                   levels = c("High-High", "Low-Low", "High-Low", 
                             "Low-High", "Not significant")))
)

# =============================================================================
# 5. GLOBAL MORAN'S I (POOLED)
# =============================================================================

message("\n=== Global Moran's I (CFR, 2007–2018 pooled) ===\n")

lw <- nb2listw(nb_2018, style = "W", zero.policy = TRUE)

# Prepare CFR vector (NA → 0 for zero-case municipalities)
cfr_for_moran <- ifelse(is.na(shp_data$cfr), 0, shp_data$cfr)

# Global Moran's I test
gm <- moran.test(cfr_for_moran, lw, zero.policy = TRUE, alternative = "greater")

message("  Moran's I = ", round(gm$estimate["Moran I statistic"], 4))
message("  p-value = ", format(gm$p.value, scientific = TRUE, digits = 3))
message("  Interpretation: ",
        ifelse(gm$p.value < 0.05,
               "SIGNIFICANT spatial clustering (HH & LL clusters exist)",
               "NO significant clustering (CFR random/independent)"))

# =============================================================================
# 6. IDENTIFY AND RANK HOTSPOTS, COLDSPOTS, OUTLIERS
# =============================================================================

message("\n=== CFR Hotspots (High-High, p < 0.05) ===\n")

hh_munis <- lisa_cfr |>
  filter(cluster_type == "High-High") |>
  arrange(desc(local_i)) |>
  select(muni_resid, local_i, p_value, cfr_pct, total_cases)

if (nrow(hh_munis) > 0) {
  message("Municipalities with significant High-High CFR clusters:\n")
  print(hh_munis)
} else {
  message("No significant High-High clusters detected.\n")
}

message("\n=== CFR Coldspots (Low-Low, p < 0.05) ===\n")

ll_munis <- lisa_cfr |>
  filter(cluster_type == "Low-Low") |>
  arrange(local_i) |>
  select(muni_resid, local_i, p_value, cfr_pct, total_cases)

if (nrow(ll_munis) > 0) {
  message("Municipalities with significant Low-Low CFR clusters:\n")
  print(ll_munis)
} else {
  message("No significant Low-Low clusters detected.\n")
}

message("\n=== CFR Spatial Outliers ===\n")

outliers_munis <- lisa_cfr |>
  filter(cluster_type %in% c("High-Low", "Low-High")) |>
  arrange(cluster_type, desc(local_i)) |>
  select(muni_resid, cluster_type, local_i, p_value, cfr_pct, total_cases)

if (nrow(outliers_munis) > 0) {
  message("Municipalities with significant spatial outlier status:\n")
  print(outliers_munis)
} else {
  message("No significant spatial outliers detected.\n")
}

# =============================================================================
# 7. ATTACH LISA RESULTS TO SHAPEFILE FOR MAPPING
# =============================================================================

map_cfr <- muni_2018 |>
  left_join(lisa_cfr, by = c("code6" = "muni_resid"))

# =============================================================================
# 8. SAVE
# =============================================================================

saveRDS(lisa_cfr,  file.path(proc_dir, "lisa_cfr_pooled.rds"))
saveRDS(map_cfr,   file.path(proc_dir, "map_lisa_cfr_pooled.rds"))

message("\nSaved:")
message("  lisa_cfr_pooled.rds      (LISA results for pooled 2007–2018)")
message("  map_lisa_cfr_pooled.rds  (Shapefile with LISA results for mapping)")

message("\nNotes:")
message("  • This LISA analysis uses POOLED 2007–2018 CFR data,")
message("    matching your INLA RQ2 approach (Models A & B)")
message("  • Results can be directly compared with INLA spatial effects (φ, u_i)")
message("  • There is NO temporal variation in this analysis")
message("    (unlike RQ1 which has RW1 temporal trend)")
message("\nNext: run 10_lisa_cfr_maps_pooled.R for figures")
