# =============================================================================
# LISA & Moran's I Maps and Figures — VL Case Fatality Rate (POOLED)
# CORRECTED: 2007–2018 pooled (matching INLA RQ2 structure)
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
fig_dir  <- "D:/AIMS/Research phase/R/figures/lisa_cfr"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD RESULTS
# =============================================================================

lisa_cfr <- readRDS(file.path(proc_dir, "lisa_cfr_pooled.rds"))
map_cfr  <- readRDS(file.path(proc_dir, "map_lisa_cfr_pooled.rds"))

message("Loaded LISA and mapping data for CFR analysis (pooled 2007–2018)")

# =============================================================================
# 2. THEME & COLORS
# =============================================================================

theme_map <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank()
    )
}

# LISA cluster colour scheme
lisa_colors <- c(
  "High-High"       = "#d73027",      # red (hotspot)
  "Low-Low"         = "#4575b4",      # blue (coldspot)
  "High-Low"        = "#fee090",      # light yellow (outlier)
  "Low-High"        = "#91bfdb",      # light blue (outlier)
  "Not significant" = "#f0f0f0"       # light grey
)

# =============================================================================
# 3. LISA CLUSTER MAP (2007–2018 POOLED)
# =============================================================================

message("\n=== Creating LISA cluster map ===")

p_lisa <- ggplot(map_cfr) +
  geom_sf(aes(fill = cluster_type), color = NA) +
  scale_fill_manual(
    values = lisa_colors,
    breaks = names(lisa_colors),
    na.value = "#f0f0f0"
  ) +
  labs(
    title = "Case Fatality Rate: LISA Cluster Map",
    subtitle = "2007–2018 (pooled)"
  ) +
  theme_map()

ggsave(
  filename = file.path(fig_dir, "lisa_cfr_pooled_2007_2018.png"),
  plot = p_lisa,
  width = 10, height = 8, dpi = 300
)

message("  Saved: lisa_cfr_pooled_2007_2018.png")

# =============================================================================
# 4. MORAN SCATTER PLOT
# =============================================================================

message("\n=== Creating Moran scatter plot ===")

scatter_data <- lisa_cfr |>
  mutate(cfr_pct = cfr_pct)

p_scatter <- ggplot(scatter_data, 
                    aes(x = cfr_pct, 
                        y = lag_cfr * 100,
                        color = cluster_type)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = mean(scatter_data$cfr_pct, na.rm = TRUE),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = mean(scatter_data$lag_cfr, na.rm = TRUE) * 100,
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(
    values = lisa_colors,
    breaks = names(lisa_colors)
  ) +
  labs(
    title = "Moran Scatter Plot: CFR",
    subtitle = "2007–2018 (pooled)",
    x = "Case Fatality Rate (%)",
    y = "Lagged Case Fatality Rate (%)",
    color = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "moran_scatter_cfr_pooled_2007_2018.png"),
  plot = p_scatter,
  width = 8, height = 7, dpi = 300
)

message("  Saved: moran_scatter_cfr_pooled_2007_2018.png")

# =============================================================================
# 5. SUMMARY TABLES
# =============================================================================

message("\n=== Summary Tables ===\n")

# High-High clusters (hotspots)
hh_summary <- lisa_cfr |>
  filter(cluster_type == "High-High") |>
  select(muni_resid, local_i, p_value, cfr_pct, total_cases) |>
  arrange(desc(local_i)) |>
  rename(
    Municipality = muni_resid,
    "Local I" = local_i,
    "p-value" = p_value,
    "CFR (%)" = cfr_pct,
    "Cases" = total_cases
  )

# Low-Low clusters (coldspots)
ll_summary <- lisa_cfr |>
  filter(cluster_type == "Low-Low") |>
  select(muni_resid, local_i, p_value, cfr_pct, total_cases) |>
  arrange(local_i) |>
  rename(
    Municipality = muni_resid,
    "Local I" = local_i,
    "p-value" = p_value,
    "CFR (%)" = cfr_pct,
    "Cases" = total_cases
  )

message("High-High municipalities (CFR hotspots):")
print(hh_summary)

message("\nLow-Low municipalities (CFR coldspots):")
print(ll_summary)

# =============================================================================
# 6. CLUSTER SUMMARY TABLE
# =============================================================================

message("\n=== Cluster Summary (2007–2018) ===\n")

cluster_summary <- lisa_cfr |>
  group_by(cluster_type) |>
  summarise(
    n = n(),
    mean_cfr = round(mean(cfr_pct, na.rm = TRUE), 2),
    median_cfr = round(median(cfr_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  arrange(factor(cluster_type,
                 levels = c("High-High", "Low-Low", "High-Low",
                           "Low-High", "Not significant")))

print(cluster_summary)

# =============================================================================
# 7. COMPARISON WITH INLA SPATIAL EFFECTS
# =============================================================================

message("\n=== Linking LISA to INLA RQ2 ===\n")
message("LISA hotspots (High-High municipalities) should correspond to:")
message("  • INLA Model A: municipalities with posterior CFR > national CFR + credible interval")
message("  • INLA Model B: municipalities with positive spatial random effect (u_i) + v_i")
message("\nWhen comparing figures:")
message("  • LISA HH (red) should overlap with INLA high RR/CFR regions")
message("  • Geographic alignment = spatial clustering confirmed in both methods")
message("  • Misalignment = outlier or boundary effect (investigate further)")

# =============================================================================
# 8. SAVE SESSION OUTPUT
# =============================================================================

message("\n=== Figure Summary ===")
message("Saved to: ", fig_dir)
message("  ✓ lisa_cfr_pooled_2007_2018.png (300 dpi)")
message("  ✓ moran_scatter_cfr_pooled_2007_2018.png (300 dpi)")

message("\nNotes:")
message("  • This is a POOLED analysis (single CFR estimate per municipality)")
message("  • Direct comparison with INLA RQ2 is now appropriate")
message("  • Temporal variation is NOT modeled (unlike RQ1)")
message("\nNext: Integrate figures & tables into thesis Chapter 4 (Results)")
