# =============================================================================
# Municipality Name Lookup
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Builds a lookup table: code6 -> municipality name + state abbreviation
# Sources: IBGE population file (contains NM_MUN or equivalent name column)
#          or the 2018 shapefile attribute table (NM_MUNICIP)
#
# Usage: source("muni_names.R") then left_join by code6 / muni_resid
#
# Output object: muni_lookup
#   columns: code6, muni_name, uf
# =============================================================================

library(sf)
library(dplyr)

shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]

shp_raw <- st_read(shp_path18, quiet = TRUE)

# Column names in the 2018 shapefile
# CD_GEOCMU = 7-digit code, NM_MUNICIP = name, NM_ESTADO = state name
# Some editions use NM_MUN instead of NM_MUNICIP — check both
name_col <- if ("NM_MUNICIP" %in% names(shp_raw)) "NM_MUNICIP" else
  if ("NM_MUN"     %in% names(shp_raw)) "NM_MUN"     else
    stop("Cannot find municipality name column in 2018 shapefile. ",
         "Available columns: ", paste(names(shp_raw), collapse = ", "))

state_names <- c(
  "11"="RO","12"="AC","13"="AM","14"="RR","15"="PA",
  "16"="AP","17"="TO","21"="MA","22"="PI","23"="CE",
  "24"="RN","25"="PB","26"="PE","27"="AL","28"="SE",
  "29"="BA","31"="MG","32"="ES","33"="RJ","35"="SP",
  "41"="PR","42"="SC","43"="RS","50"="MS","51"="MT",
  "52"="GO","53"="DF"
)

muni_lookup <- shp_raw |>
  st_drop_geometry() |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(
    code6      = substr(as.character(CD_GEOCMU), 1, 6),
    muni_name  = as.character(.data[[name_col]]),
    state_code = substr(code6, 1, 2),
    uf         = state_names[state_code]
  ) |>
  select(code6, muni_name, uf) |>
  distinct()

message("Municipality lookup built: ", nrow(muni_lookup), " entries")