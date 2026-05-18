# =============================================================================
# Publication-Quality Maps — VL Spatial & Spatio-Temporal Analysis
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Maps produced:
#   MAP 1: Incidence rate choropleth — 3 period panels (2001-06, 2007-12, 2013-18)
#   MAP 2: LISA cluster map          — 3 period panels
#   MAP 3: Annual HH cluster count   — time series bar chart
#   MAP 4: Smoothed incidence rate   — spotlight years 2001, 2009, 2018
#
# Packages: ggplot2, sf, dplyr, patchwork, scales, RColorBrewer
#   install.packages(c("ggplot2","patchwork","scales","RColorBrewer"))
#
# Output: PNG files at 300 DPI in outputs/maps/
#
# Geometry note: 2018 shapefile contains one invalid polygon in Pará (state 15).
#   Fix applied: sf_use_s2(FALSE) + st_make_valid() inside make_state_borders().
#
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
out_dir  <- "D:/AIMS/Research phase/R/outputs/maps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

map_p1       <- readRDS(file.path(proc_dir, "map_lisa_2001_2006.rds"))
map_p2       <- readRDS(file.path(proc_dir, "map_lisa_2007_2012.rds"))
map_p3       <- readRDS(file.path(proc_dir, "map_lisa_2013_2018.rds"))
agg_smoothed <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))
lisa_annual  <- readRDS(file.path(proc_dir, "lisa_annual.rds"))

shp_dir07 <- "data/shapefiles/2007"
shp_dir18 <- "data/shapefiles/2018"

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

message("Data loaded.")

# =============================================================================
# 2. SHARED HELPERS
# =============================================================================

# State border polygons (for reference overlay on maps)
# sf_use_s2(FALSE) required: 2018 shapefile has one invalid polygon in Pará
make_state_borders <- function(shp) {
  sf_use_s2(FALSE)
  result <- shp |>
    st_make_valid() |>
    mutate(state_code = substr(code6, 1, 2)) |>
    group_by(state_code) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
  sf_use_s2(TRUE)
  result
}

message("Building state borders...")
state_borders_07 <- make_state_borders(muni_2007)
state_borders_18 <- make_state_borders(muni_2018)
message("  ", nrow(state_borders_07), " states (2007), ",
        nrow(state_borders_18), " states (2018)")

# Clean map theme for publication
theme_vl_map <- function(...) {
  theme_void(...) +
    theme(
      plot.title        = element_text(size = 11, face = "bold",
                                       hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle     = element_text(size = 9, hjust = 0.5,
                                       color = "grey40", margin = margin(b = 6)),
      legend.position   = "bottom",
      legend.title      = element_text(size = 8, face = "bold"),
      legend.text       = element_text(size = 7),
      legend.key.width  = unit(1.2, "cm"),
      legend.key.height = unit(0.35, "cm"),
      plot.margin       = margin(4, 4, 4, 4)
    )
}

choro_palette <- c("grey92", brewer.pal(9, "YlOrRd"))

# =============================================================================
# 3. MAP 1 — INCIDENCE RATE CHOROPLETH (3 PERIODS)
# =============================================================================

message("Building Map 1: Incidence choropleth...")

make_period_rate_sf <- function(years, shp) {
  rates <- agg_smoothed |>
    filter(year %in% years) |>
    group_by(muni_resid) |>
    summarise(
      total_cases = sum(n_cases,   na.rm = TRUE),
      total_pop   = sum(POPULACAO, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      period_rate  = ifelse(total_pop > 0,
                            total_cases / total_pop * 100000, NA_real_),
      rate_display = pmin(period_rate,
                          quantile(period_rate, 0.99, na.rm = TRUE))
    )
  shp |> left_join(rates, by = c("code6" = "muni_resid"))
}

choro_p1 <- make_period_rate_sf(2001:2006, muni_2007)
choro_p2 <- make_period_rate_sf(2007:2012, muni_2018)
choro_p3 <- make_period_rate_sf(2013:2018, muni_2018)

rate_max <- max(
  quantile(choro_p1$period_rate, 0.99, na.rm = TRUE),
  quantile(choro_p2$period_rate, 0.99, na.rm = TRUE),
  quantile(choro_p3$period_rate, 0.99, na.rm = TRUE)
)

make_choro_panel <- function(data, state_borders, title, subtitle,
                             rate_max, palette) {
  ggplot(data) +
    geom_sf(aes(fill = rate_display), colour = NA, linewidth = 0) +
    geom_sf(data = state_borders, fill = NA,
            colour = "white", linewidth = 0.25) +
    scale_fill_gradientn(
      colours  = palette,
      na.value = "grey85",
      limits   = c(0, rate_max),
      name     = "Incidence rate\n(per 100,000)",
      labels   = label_number(accuracy = 1),
      guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                 barwidth = 8, barheight = 0.5)
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_vl_map()
}

p_choro1 <- make_choro_panel(choro_p1, state_borders_07,
                             "2001–2006", "Baseline period",
                             rate_max, choro_palette)
p_choro2 <- make_choro_panel(choro_p2, state_borders_18,
                             "2007–2012", "Expansion phase",
                             rate_max, choro_palette)
p_choro3 <- make_choro_panel(choro_p3, state_borders_18,
                             "2013–2018", "Recent period",
                             rate_max, choro_palette)

map1 <- (p_choro1 | p_choro2 | p_choro3) +
  plot_annotation(
    title    = "Visceral Leishmaniasis Incidence Rate by Municipality, Brazil 2001–2018",
    subtitle = "Global empirical Bayesian smoothed rates, capped at 99th percentile",
    caption  = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas.",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
      plot.caption  = element_text(size = 7,  hjust = 1,   color = "grey50")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "map1_incidence_choropleth.png"),
       map1, width = 14, height = 7, dpi = 300, bg = "white")
message("  Map 1 saved.")

# =============================================================================
# 4. MAP 2 — LISA CLUSTER MAP (3 PERIODS)
# =============================================================================

message("Building Map 2: LISA cluster maps...")

lisa_colours <- c(
  "High-High"       = "#d7191c",
  "Low-Low"         = "#2c7bb6",
  "High-Low"        = "#fdae61",
  "Low-High"        = "#abd9e9",
  "Not significant" = "grey90"
)

lisa_labels <- c(
  "High-High"       = "High-High (hotspot)",
  "Low-Low"         = "Low-Low (coldspot)",
  "High-Low"        = "High-Low (outlier)",
  "Low-High"        = "Low-High (outlier)",
  "Not significant" = "Not significant"
)

make_lisa_panel <- function(map_sf, state_borders, title, subtitle) {
  map_sf <- map_sf |>
    mutate(cluster_type = factor(
      coalesce(cluster_type, "Not significant"),
      levels = names(lisa_colours)
    ))
  ggplot(map_sf) +
    geom_sf(aes(fill = cluster_type), colour = NA, linewidth = 0) +
    geom_sf(data = state_borders, fill = NA,
            colour = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = lisa_colours,
      labels = lisa_labels,
      name   = "LISA cluster type",
      drop   = FALSE,
      guide  = guide_legend(title.position = "top", title.hjust = 0.5,
                            nrow = 2,
                            override.aes = list(colour = "grey70",
                                                linewidth = 0.2))
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_vl_map()
}

p_lisa1 <- make_lisa_panel(map_p1, state_borders_07,
                           "2001–2006", "Baseline period")
p_lisa2 <- make_lisa_panel(map_p2, state_borders_18,
                           "2007–2012", "Expansion phase")
p_lisa3 <- make_lisa_panel(map_p3, state_borders_18,
                           "2013–2018", "Recent period")

map2 <- (p_lisa1 | p_lisa2 | p_lisa3) +
  plot_annotation(
    title    = "LISA Spatial Cluster Map — Visceral Leishmaniasis, Brazil 2001–2018",
    subtitle = paste0("Local Moran's I, p < 0.05 (999 permutations). ",
                      "Global Moran's I: 0.246 (2001-06), 0.230 (2007-12), 0.090 (2013-18)"),
    caption  = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas.",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
      plot.caption  = element_text(size = 7,  hjust = 1,   color = "grey50")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "map2_LISA_clusters.png"),
       map2, width = 14, height = 7, dpi = 300, bg = "white")
message("  Map 2 saved.")

# =============================================================================
# 5. MAP 3 — ANNUAL HH CLUSTER COUNT TIME SERIES
# =============================================================================

message("Building Map 3: Annual HH cluster count...")

hh_annual <- lisa_annual |>
  group_by(year) |>
  summarise(n_HH = sum(cluster_type == "High-High"), .groups = "drop")

map3 <- ggplot(hh_annual, aes(x = year, y = n_HH)) +
  geom_col(fill = "#d7191c", alpha = 0.85, width = 0.75) +
  geom_line(colour = "#8b0000", linewidth = 0.8) +
  geom_point(colour = "#8b0000", size = 2.5) +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Number of High-High Spatial Clusters per Year",
    subtitle = paste0("Municipalities with significantly elevated VL incidence\n",
                      "surrounded by high-incidence neighbours (Local Moran's I, p < 0.05)"),
    x = "Year", y = "Number of High-High municipalities",
    caption = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
    plot.caption       = element_text(size = 7, color = "grey50", hjust = 1),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "map3_annual_HH_clusters.png"),
       map3, width = 10, height = 5, dpi = 300, bg = "white")
message("  Map 3 saved.")

# =============================================================================
# 6. MAP 4 — SPOTLIGHT MAPS: SELECTED YEARS 2001, 2009, 2018
# =============================================================================

message("Building Map 4: Spotlight year maps...")

spotlight_max <- quantile(
  agg_smoothed$smoothed_rate[agg_smoothed$smoothed_rate > 0],
  0.99, na.rm = TRUE
)

make_spotlight <- function(yr, shp, state_borders, rate_max) {
  yr_data <- agg_smoothed |>
    filter(year == yr) |>
    select(muni_resid, smoothed_rate) |>
    mutate(rate_display = pmin(smoothed_rate, rate_max))
  sf_data <- shp |>
    left_join(yr_data, by = c("code6" = "muni_resid"))
  ggplot(sf_data) +
    geom_sf(aes(fill = rate_display), colour = NA, linewidth = 0) +
    geom_sf(data = state_borders, fill = NA,
            colour = "white", linewidth = 0.25) +
    scale_fill_gradientn(
      colours  = c("grey92", brewer.pal(9, "YlOrRd")),
      na.value = "grey85",
      limits   = c(0, rate_max),
      name     = "Rate per 100,000",
      labels   = label_number(accuracy = 1),
      guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                 barwidth = 6, barheight = 0.4)
    ) +
    labs(title = as.character(yr)) +
    theme_vl_map()
}

sp_2001 <- make_spotlight(2001, muni_2007, state_borders_07, spotlight_max)
sp_2009 <- make_spotlight(2009, muni_2007, state_borders_07, spotlight_max)
sp_2018 <- make_spotlight(2018, muni_2018, state_borders_18, spotlight_max)

map4 <- (sp_2001 | sp_2009 | sp_2018) +
  plot_annotation(
    title    = "Smoothed VL Incidence Rate: Selected Years 2001, 2009, 2018",
    subtitle = "Global empirical Bayesian smoothed rates per 100,000 inhabitants",
    caption  = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas.",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
      plot.caption  = element_text(size = 7,  hjust = 1,   color = "grey50")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "map4_spotlight_years.png"),
       map4, width = 14, height = 7, dpi = 300, bg = "white")
message("  Map 4 saved.")

# =============================================================================
# 7. SUMMARY
# =============================================================================

message("\n=== All maps saved to: ", out_dir, " ===")
message("  map1_incidence_choropleth.png  — 3-period incidence rate panels")
message("  map2_LISA_clusters.png         — 3-period LISA cluster panels")
message("  map3_annual_HH_clusters.png    — annual HH count bar chart")
message("  map4_spotlight_years.png       — smoothed rates for 2001, 2009, 2018")
message("\nAll maps: 300 DPI, suitable for thesis submission.")