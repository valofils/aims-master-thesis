# =============================================================================
# Local Indicators of Spatial Association (LISA) — VL Incidence Rates
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Cluster types:
#   High-High (HH): hotspot — high rate surrounded by high-rate neighbours
#   Low-Low   (LL): coldspot — low rate surrounded by low-rate neighbours
#   High-Low  (HL): spatial outlier — isolated high rate
#   Low-High  (LH): spatial outlier — isolated low rate
#   Not significant: p >= 0.05
#
# Three analysis periods (following de Melo et al. 2023, Bruhn et al. 2024):
#   Period 1: 2001–2006  baseline
#   Period 2: 2007–2012  expansion
#   Period 3: 2013–2018  recent
#
# spdep version: 1.4.2
#   localmoran() no longer accepts nsim — use localmoran_perm() for
#   permutation-based p-values.
#
# =============================================================================

library(sf)
library(dplyr)
library(spdep)
library(tidyr)

proc_dir  <- "D:/AIMS/Research phase/R/data/processed"
shp_dir07 <- "data/shapefiles/2007"
shp_dir18 <- "data/shapefiles/2018"

set.seed(2024)

# =============================================================================
# 1. LOAD DATA AND SHAPEFILES
# =============================================================================

agg_smoothed <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))
nb_2007      <- readRDS(file.path(proc_dir, "nb_2007.rds"))
nb_2018      <- readRDS(file.path(proc_dir, "nb_2018.rds"))

message("Aggregated smoothed data: ", nrow(agg_smoothed), " rows")

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
# 2. CORE LISA FUNCTION — spdep 1.4.2 compatible
# =============================================================================

run_lisa <- function(rate_vec, nb, shapefile, label = "") {
  
  n <- length(rate_vec)
  stopifnot(length(nb) == n)
  
  # Replace NA with 0 — missing-population municipalities cannot form clusters
  rate_clean <- ifelse(is.na(rate_vec), 0, rate_vec)
  
  # Row-standardised spatial weights
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # localmoran_perm: permutation-based p-values (spdep >= 1.2)
  # Falls back to localmoran (analytical) if localmoran_perm unavailable
  local_moran <- tryCatch(
    localmoran_perm(
      x           = rate_clean,
      listw       = lw,
      zero.policy = TRUE,
      alternative = "two.sided",
      nsim        = 999,
      iseed       = 2024
    ),
    error = function(e) {
      message("  Falling back to localmoran() (analytical): ", e$message)
      localmoran(
        x           = rate_clean,
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
  x_mean <- mean(rate_clean, na.rm = TRUE)
  lag_x  <- lag.listw(lw, rate_clean, zero.policy = TRUE)
  
  cluster_type <- case_when(
    p_value >= 0.05                        ~ "Not significant",
    rate_clean >= x_mean & lag_x >= x_mean ~ "High-High",
    rate_clean <  x_mean & lag_x <  x_mean ~ "Low-Low",
    rate_clean >= x_mean & lag_x <  x_mean ~ "High-Low",
    rate_clean <  x_mean & lag_x >= x_mean ~ "Low-High",
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
    lag_rate     = lag_x,
    cluster_type = cluster_type,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 3. PERIOD-LEVEL LISA
# =============================================================================

message("\n=== Period-level LISA ===")

periods <- list(
  list(label = "2001-2006", years = 2001:2006, shp = muni_2007, nb = nb_2007),
  list(label = "2007-2012", years = 2007:2012, shp = muni_2018, nb = nb_2018),
  list(label = "2013-2018", years = 2013:2018, shp = muni_2018, nb = nb_2018)
)

lisa_period_list <- vector("list", 3)

for (p in seq_along(periods)) {
  
  pd    <- periods[[p]]
  label <- pd$label
  shp   <- pd$shp
  nb    <- pd$nb
  
  message("\nPeriod: ", label)
  
  # Aggregate cases and population across period years
  period_agg <- agg_smoothed |>
    filter(year %in% pd$years) |>
    group_by(muni_resid) |>
    summarise(
      total_cases = sum(n_cases,   na.rm = TRUE),
      total_pop   = sum(POPULACAO, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      period_rate = ifelse(total_pop > 0,
                           total_cases / total_pop * 100000,
                           NA_real_)
    )
  
  # Align to shapefile row order
  shp_data <- shp |>
    st_drop_geometry() |>
    left_join(period_agg, by = c("code6" = "muni_resid"))
  
  # Run LISA
  lisa_res <- run_lisa(shp_data$period_rate, nb, shp, label = label)
  lisa_res$period      <- label
  lisa_res$total_cases <- shp_data$total_cases
  lisa_res$period_rate <- shp_data$period_rate
  
  lisa_period_list[[p]] <- lisa_res
}

lisa_periods <- bind_rows(lisa_period_list)

message("\n--- Period LISA cluster counts ---")
print(
  lisa_periods |>
    group_by(period, cluster_type) |>
    summarise(n = n(), .groups = "drop") |>
    pivot_wider(names_from = cluster_type, values_from = n,
                values_fill = 0)
)

# =============================================================================
# 4. ANNUAL LISA
# =============================================================================

message("\n=== Annual LISA ===")

lisa_annual_list <- vector("list", 18)

for (yr in 2001:2018) {
  
  shp <- if (yr <= 2009) muni_2007 else muni_2018
  nb  <- if (yr <= 2009) nb_2007   else nb_2018
  
  shp_data <- shp |>
    st_drop_geometry() |>
    left_join(
      agg_smoothed |>
        filter(year == yr) |>
        select(muni_resid, smoothed_rate, n_cases),
      by = c("code6" = "muni_resid")
    )
  
  lisa_res <- run_lisa(shp_data$smoothed_rate, nb, shp,
                       label = as.character(yr))
  lisa_res$year    <- yr
  lisa_res$n_cases <- shp_data$n_cases
  
  lisa_annual_list[[yr - 2000]] <- lisa_res
}

lisa_annual <- bind_rows(lisa_annual_list)

message("\n--- Annual HH cluster counts ---")
print(
  lisa_annual |>
    filter(cluster_type == "High-High") |>
    group_by(year) |>
    summarise(n_HH = n(), .groups = "drop"),
  n = 18
)

# =============================================================================
# 5. GLOBAL MORAN'S I PER PERIOD
# =============================================================================

message("\n=== Global Moran's I per period ===")

for (p in seq_along(periods)) {
  
  pd    <- periods[[p]]
  shp   <- pd$shp
  nb    <- pd$nb
  label <- pd$label
  
  period_agg <- agg_smoothed |>
    filter(year %in% pd$years) |>
    group_by(muni_resid) |>
    summarise(total_cases = sum(n_cases, na.rm = TRUE),
              total_pop   = sum(POPULACAO, na.rm = TRUE),
              .groups = "drop") |>
    mutate(period_rate = ifelse(total_pop > 0,
                                total_cases / total_pop * 100000, 0))
  
  shp_data <- shp |>
    st_drop_geometry() |>
    left_join(period_agg, by = c("code6" = "muni_resid")) |>
    mutate(period_rate = coalesce(period_rate, 0))
  
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  gm <- moran.test(shp_data$period_rate, lw,
                   zero.policy = TRUE, alternative = "greater")
  
  message("  Period ", label,
          ":  Moran's I = ",
          round(gm$estimate["Moran I statistic"], 4),
          "  p = ", format(gm$p.value, scientific = TRUE, digits = 3))
}

# =============================================================================
# 6. ATTACH LISA RESULTS TO SHAPEFILES FOR MAPPING
# =============================================================================

map_p1 <- muni_2007 |>
  left_join(lisa_periods |> filter(period == "2001-2006"),
            by = c("code6" = "muni_resid"))

map_p2 <- muni_2018 |>
  left_join(lisa_periods |> filter(period == "2007-2012"),
            by = c("code6" = "muni_resid"))

map_p3 <- muni_2018 |>
  left_join(lisa_periods |> filter(period == "2013-2018"),
            by = c("code6" = "muni_resid"))

# =============================================================================
# 7. SAVE
# =============================================================================

saveRDS(lisa_periods, file.path(proc_dir, "lisa_periods.rds"))
saveRDS(lisa_annual,  file.path(proc_dir, "lisa_annual.rds"))
saveRDS(map_p1,       file.path(proc_dir, "map_lisa_2001_2006.rds"))
saveRDS(map_p2,       file.path(proc_dir, "map_lisa_2007_2012.rds"))
saveRDS(map_p3,       file.path(proc_dir, "map_lisa_2013_2018.rds"))

message("\nSaved:")
message("  lisa_periods.rds")
message("  lisa_annual.rds")
message("  map_lisa_2001_2006.rds  |  map_lisa_2007_2012.rds  |  map_lisa_2013_2018.rds")
message("\nNext: run 06_maps.R")