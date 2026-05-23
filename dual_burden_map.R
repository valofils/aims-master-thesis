# =============================================================================
# dual_burden_map.R
# Dual-Burden Municipalities: High VL Incidence and High Case Fatality
#
# Inputs (from proc_dir):
#   inla_rq1_spatial.rds   — spatial_combined with p_exceed (P(RR > 1))
#   inla_rq2_muni_cfr.rds  — post_cfr with p_exceed_cfr (P(CFR > national))
#
# Output:
#   rq2_dual_burden_map.png
# =============================================================================

library(sf)
library(ggplot2)
library(dplyr)
library(scales)

# ---- Paths ------------------------------------------------------------------

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
fig_dir  <- "D:/AIMS/Research phase/R/figures/inla_rq2"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load shapefile ---------------------------------------------------------

shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]
muni_2018  <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

# ---- State borders ----------------------------------------------------------

sf_use_s2(FALSE)
state_borders_18 <- muni_2018 |>
  st_make_valid() |>
  mutate(state_code = substr(code6, 1, 2)) |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")
sf_use_s2(TRUE)

# ---- Load INLA results ------------------------------------------------------

spatial_rq1 <- readRDS(file.path(proc_dir, "inla_rq1_spatial.rds"))
post_cfr    <- readRDS(file.path(proc_dir, "inla_rq2_muni_cfr.rds"))

# ---- Thresholds -------------------------------------------------------------

thresh_inc <- 0.95   # P(spatial RR > 1) > 0.95  → high incidence
thresh_cfr <- 0.95   # P(CFR > national CFR) > 0.95 → high CFR

# ---- Classify municipalities ------------------------------------------------

high_inc <- spatial_rq1 |>
  filter(p_exceed > thresh_inc) |>
  pull(muni_resid)

high_cfr_munis <- post_cfr |>
  filter(p_exceed_cfr > thresh_cfr) |>
  pull(muni_resid)

n_high_inc  <- length(high_inc)
n_high_cfr  <- length(high_cfr_munis)
n_dual      <- length(intersect(high_inc, high_cfr_munis))

message("High incidence municipalities (P(RR>1)>0.95): ", n_high_inc)
message("High CFR municipalities (P(CFR>nat)>0.95):    ", n_high_cfr)
message("Dual burden municipalities:                    ", n_dual)

# ---- Build national CFR threshold label ------------------------------------
# The CFR threshold shown in the subtitle is the national posterior CFR
nat_cfr_pct <- round(
  plogis(post_cfr$cfr_mean[1]) * 100,   # fallback — use the value from script
  1
)
# Retrieve it properly from model intercept if available; otherwise hard-code
# from known result: national CFR ≈ 6.7%
nat_cfr_label <- "6.7%"

# ---- Join categories to shapefile ------------------------------------------

dual_sf <- muni_2018 |>
  mutate(
    category = case_when(
      code6 %in% intersect(high_inc, high_cfr_munis) ~ "Dual burden",
      code6 %in% high_inc                             ~ "High incidence only",
      code6 %in% high_cfr_munis                       ~ "High CFR only",
      TRUE                                             ~ "Neither"
    ),
    category = factor(category,
                      levels = c("Dual burden",
                                 "High incidence only",
                                 "High CFR only",
                                 "Neither"))
  )

# ---- Colour palette ---------------------------------------------------------

cat_colours <- c(
  "Dual burden"         = "#6A0DAD",   # purple
  "High incidence only" = "#CC0000",   # red
  "High CFR only"       = "#E6A817",   # orange/amber
  "Neither"             = "#D3D3D3"    # light grey
)

# ---- Plot -------------------------------------------------------------------

p_dual <- ggplot(dual_sf) +
  geom_sf(aes(fill = category), colour = NA, linewidth = 0) +
  geom_sf(data = state_borders_18, fill = NA,
          colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values   = cat_colours,
    name     = NULL,
    drop     = FALSE,
    guide    = guide_legend(
      nrow             = 2,
      byrow            = TRUE,
      title.position   = "top",
      label.theme      = element_text(size = 9),
      keywidth         = unit(0.9, "cm"),
      keyheight        = unit(0.5, "cm")
    )
  ) +
  labs(
    caption  = "Source: SINAN/IBGE. Agusto Lab, University of Kansas."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.caption     = element_text(size = 7,  hjust = 1, colour = "grey50"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 10),
    plot.margin      = margin(10, 10, 10, 10)
  )

# ---- Save -------------------------------------------------------------------

ggsave(
  file.path(fig_dir, "rq2_dual_burden_map.png"),
  p_dual, width = 8, height = 9, dpi = 300, bg = "white"
)
message("Saved: rq2_dual_burden_map.png")
