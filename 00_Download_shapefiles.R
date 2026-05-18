# =============================================================================
# IBGE Municipal Shapefiles — 2007 and 2018 editions
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Confirmed from iterative inspection:
#
#   2007  file : 55mu2500gsr.shp   |  code col : GEOCODIG_M  (7 digits)
#   2018  file : BRMUE250GC_SIR.shp | code col : CD_GEOCMU   (7 digits)
#
#   Both files contain two identical non-municipality water-body rows:
#     4300001  Lagoa Mirim     (Rio Grande do Sul)
#     4300002  Lagoa dos Patos (Rio Grande do Sul)
#   These appear with NA names in 2007. Remove by exact 7-digit code.
#
#   Join key: first 6 digits of 7-digit code
#             matches V9 / V19 in SINAN and MUNIC_RES in population data
#
# =============================================================================

library(sf)
library(dplyr)

options(timeout = 300)

dir.create("data/shapefiles/2007", recursive = TRUE, showWarnings = FALSE)
dir.create("data/shapefiles/2018", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. DOWNLOAD
# =============================================================================

url_2007 <- paste0(
  "https://geoftp.ibge.gov.br/organizacao_do_territorio/",
  "malhas_territoriais/malhas_municipais/municipio_2007/",
  "escala_2500mil/proj_geografica_sirgas2000/brasil/",
  "55mu2500gsr.zip"
)

url_2018 <- paste0(
  "https://geoftp.ibge.gov.br/organizacao_do_territorio/",
  "malhas_territoriais/malhas_municipais/municipio_2018/",
  "Brasil/BR/",
  "br_municipios.zip"
)

zip_2007 <- "data/shapefiles/2007/55mu2500gsr.zip"
zip_2018 <- "data/shapefiles/2018/br_municipios.zip"

if (!file.exists(zip_2007)) {
  message("Downloading 2007 shapefile (8.5 MB) ...")
  download.file(url_2007, destfile = zip_2007, mode = "wb", method = "libcurl")
  message("2007 download complete.")
} else {
  message("2007 zip already present — skipping.")
}

if (!file.exists(zip_2018)) {
  message("Downloading 2018 shapefile (122 MB) ...")
  download.file(url_2018, destfile = zip_2018, mode = "wb", method = "libcurl")
  message("2018 download complete.")
} else {
  message("2018 zip already present — skipping.")
}

# =============================================================================
# 2. UNZIP AND READ
# =============================================================================

unzip(zip_2007, exdir = "data/shapefiles/2007", overwrite = TRUE)
unzip(zip_2018, exdir = "data/shapefiles/2018", overwrite = TRUE)

shp_2007_path <- list.files("data/shapefiles/2007", pattern = "\\.shp$",
                            full.names = TRUE, recursive = TRUE)
shp_2018_path <- list.files("data/shapefiles/2018", pattern = "\\.shp$",
                            full.names = TRUE, recursive = TRUE)

muni_2007_raw <- st_read(shp_2007_path[1], quiet = TRUE)
muni_2018_raw <- st_read(shp_2018_path[1], quiet = TRUE)

message("2007 raw: ", nrow(muni_2007_raw), " rows")  # 5566
message("2018 raw: ", nrow(muni_2018_raw), " rows")  # 5572

# =============================================================================
# 3. CLEAN BOTH SHAPEFILES
# =============================================================================
# Remove Lagoa Mirim and Lagoa dos Patos — present in both files, same codes.
# Then take first 6 digits of the 7-digit IBGE code as the join key.

WATER_BODIES <- c("4300001", "4300002")

muni_2007 <- muni_2007_raw |>
  filter(!as.character(GEOCODIG_M) %in% WATER_BODIES) |>
  mutate(code6 = substr(as.character(GEOCODIG_M), 1, 6))

muni_2018 <- muni_2018_raw |>
  filter(!as.character(CD_GEOCMU) %in% WATER_BODIES) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

message("2007 after cleaning: ", nrow(muni_2007), " municipalities")  # 5564
message("2007 duplicates:     ", sum(duplicated(muni_2007$code6)))     # 0

message("2018 after cleaning: ", nrow(muni_2018), " municipalities")  # 5570
message("2018 duplicates:     ", sum(duplicated(muni_2018$code6)))     # 0

# =============================================================================
# 4. UNIQUENESS CHECK
# =============================================================================

stopifnot(
  "Duplicate code6 in 2007 shapefile" = !any(duplicated(muni_2007$code6)),
  "Duplicate code6 in 2018 shapefile" = !any(duplicated(muni_2018$code6))
)
message("\nAll code6 keys unique — shapefiles ready to join.")

# =============================================================================
# 5. QUICK PLOT
# =============================================================================

par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
plot(st_geometry(muni_2007), col = "grey90", border = "white",
     lwd = 0.1, main = paste0("2007  (n = ", nrow(muni_2007), ")"))
plot(st_geometry(muni_2018), col = "grey90", border = "white",
     lwd = 0.1, main = paste0("2018  (n = ", nrow(muni_2018), ")"))
par(mfrow = c(1, 1))

# =============================================================================
# 6. COVERAGE CHECK (run after loading population data)
# =============================================================================
# pop$code6 <- substr(as.character(pop$MUNIC_RES), 1, 6)
#
# miss_07 <- setdiff(unique(pop$code6[pop$ANO <= 2009]), muni_2007$code6)
# miss_18 <- setdiff(unique(pop$code6[pop$ANO >= 2010]), muni_2018$code6)
# message("Pop codes not in 2007 shapefile: ", length(miss_07))
# message("Pop codes not in 2018 shapefile: ", length(miss_18))

# =============================================================================
# 7. USAGE
# =============================================================================
# muni_2007  ->  join to data where year %in% 2001:2009
# muni_2018  ->  join to data where year %in% 2010:2018
#
# Join key: code6  <->  V9 / V19 (SINAN)  <->  substr(MUNIC_RES, 1, 6) (pop)
#
# Example:
#   map_early <- left_join(muni_2007,
#                          cases_agg |> filter(year <= 2009),
#                          by = c("code6" = "muni_code"))
#
#   map_late  <- left_join(muni_2018,
#                          cases_agg |> filter(year >= 2010),
#                          by = c("code6" = "muni_code"))

message("Done. Objects in memory: muni_2007, muni_2018")