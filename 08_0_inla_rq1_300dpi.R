# =============================================================================
# INLA — RQ1: Spatio-Temporal Model of VL Incidence
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Model:
#   Y_it ~ Poisson(E_it * theta_it)
#   log(theta_it) = alpha + u_i + v_i + gamma_t
#
#   where:
#     Y_it       = observed confirmed VL cases, municipality i, year t
#     E_it       = expected cases (population offset using national rate)
#     alpha      = intercept (overall log relative risk)
#     u_i + v_i  = BYM2 spatial random effect (Riebler et al. 2016)
#                  u_i = spatially structured (ICAR, Besag)
#                  v_i = unstructured (iid)
#     gamma_t    = RW1 temporal random effect over years 2001–2018
#
# Spatial graph: nb_2018.graph (5,570 municipalities, 2018 shapefile)
#   Used consistently across all 18 years as the single spatial domain.
#
# Outputs saved to processed/:
#   inla_rq1_model.rds       — full INLA model object
#   inla_rq1_results.rds     — tidy results data frame (municipality-year RR)
#   inla_rq1_temporal.rds    — posterior temporal trend (gamma_t)
#   inla_rq1_spatial.rds     — posterior spatial effects per municipality
#
# Figures saved to figures/inla_rq1/:
#   rq1_temporal_trend.png   — posterior RW1 trend with 95% credible interval
#   rq1_rr_map_period.png    — posterior mean RR map (3 period averages)
#   rq1_exceedance_map.png   — P(RR > 1) exceedance probability map
#
# References:
#   Rue, Martino & Chopin (2009) JRSS-B — R-INLA
#   Riebler et al. (2016) Statistical Methods in Medical Research — BYM2
#   Moraga (2019) Geospatial Health Data — BYM implementation in INLA
#   Clayton & Kaldor (1987) Biometrics — empirical Bayes smoothing
#
# =============================================================================

library(INLA)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
fig_dir  <- "D:/AIMS/Research phase/R/figures/inla_rq1"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("Loading data...")

agg <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))

message("  Rows: ", nrow(agg))
message("  Years: ", min(agg$year), "–", max(agg$year))
message("  Municipalities: ", n_distinct(agg$muni_resid))

# Load shapefile for mapping outputs
shp_dir18   <- "data/shapefiles/2018"
shp_path18  <- list.files(shp_dir18, pattern = "\\.shp$",
                          full.names = TRUE, recursive = TRUE)[1]
muni_2018   <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

graph_path <- file.path(proc_dir, "nb_2018.graph")
stopifnot("nb_2018.graph not found — run nb2INLA() first" = file.exists(graph_path))

# =============================================================================
# 2. PREPARE INLA DATA FRAME
# =============================================================================

message("Preparing INLA data frame...")

# National rate (per person) used to compute expected counts E_it
# E_it = pop_it * national_rate  =>  if theta=1 everywhere, sum(Y) = sum(E)
national_rate <- sum(agg$n_cases, na.rm = TRUE) /
  sum(agg$POPULACAO, na.rm = TRUE)

message("  National rate (per person): ", formatC(national_rate, format = "e", digits = 4))

# Municipality index: integer 1..N aligned to the graph node order.
# The graph was built from nb_2018 which is ordered by shapefile row.
# We use the same ordering here via match().
muni_order <- muni_2018$code6   # 5570 codes in shapefile row order

# Year index: integer 1..18  (2001 = 1, ..., 2018 = 18)
inla_df <- agg |>
  mutate(
    E_it       = POPULACAO * national_rate,   # expected count
    # INLA does not allow NA in E — replace with small positive value and
    # set n_cases to NA so INLA treats these rows as missing observations
    n_cases    = ifelse(is.na(POPULACAO) | POPULACAO <= 0, NA_integer_, n_cases),
    E_it       = ifelse(is.na(E_it) | E_it <= 0, 1e-6, E_it),
    muni_idx   = match(muni_resid, muni_order),  # 1..5570
    year_idx   = year - 2000L,                   # 1..18
    muni_idx2  = muni_idx
  ) |>
  filter(!is.na(muni_idx))   # drop any codes not in 2018 shapefile

message("  Rows after index alignment: ", nrow(inla_df))
message("  Municipalities matched:     ", n_distinct(inla_df$muni_idx))
message("  Expected count range: [",
        round(min(inla_df$E_it, na.rm = TRUE), 3), ", ",
        round(max(inla_df$E_it, na.rm = TRUE), 1), "]")

# Check: rows with NA expected count
n_na_E <- sum(is.na(inla_df$E_it))
message("  Rows with NA E_it (missing population): ", n_na_E)

# =============================================================================
# 3. PRIORS
# =============================================================================
# BYM2 prior (Riebler et al. 2016):
#   - precision of spatial effect: PC prior, P(sigma > 1) = 0.01
#   - mixing parameter phi (proportion spatially structured):
#     PC prior, P(phi < 0.5) = 0.5  (weakly informative, no strong assumption)
#
# RW1 temporal prior:
#   - precision: PC prior, P(sigma > 0.5) = 0.01  (smooth trend)
#
# These are weakly informative and standard for disease mapping.
# See Moraga (2019) chapter 6 and Riebler et al. (2016).

prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

prior_rw1 <- list(
  prec = list(prior = "pc.prec", param = c(0.5, 0.01))
)

# =============================================================================
# 4. INLA FORMULA
# =============================================================================

formula_rq1 <- n_cases ~
  1 +                                           # intercept = alpha
  f(muni_idx,                                   # BYM2 spatial random effect
    model        = "bym2",
    graph        = graph_path,
    scale.model  = TRUE,
    constr       = TRUE,
    hyper        = prior_bym2) +
  f(year_idx,                                   # RW1 temporal random effect
    model        = "rw1",
    scale.model  = TRUE,
    constr       = TRUE,
    hyper        = prior_rw1)

# =============================================================================
# 5. FIT MODEL
# =============================================================================

message("\nFitting INLA Poisson BYM2 + RW1 model...")
message("  This may take 5–15 minutes depending on your machine.")

t_start <- proc.time()

model_rq1 <- inla(
  formula  = formula_rq1,
  family   = "poisson",
  data     = inla_df,
  E        = E_it,                        # expected counts offset
  control.compute = list(
    dic        = TRUE,
    waic       = TRUE,
    cpo        = TRUE,
    config     = TRUE,                    # needed for posterior sampling
    return.marginals.predictor = FALSE    # saves memory
  ),
  control.predictor = list(
    compute = TRUE,
    link    = 1
  ),
  control.inla = list(
    strategy   = "adaptive",
    int.strategy = "eb"                   # empirical Bayes integration (faster)
  ),
  verbose = FALSE
)

t_elapsed <- proc.time() - t_start
message("  Model fitted in ", round(t_elapsed["elapsed"] / 60, 1), " minutes.")

# =============================================================================
# 6. MODEL FIT SUMMARY
# =============================================================================

message("\n=== Model fit summary ===")
message("  DIC:  ", round(model_rq1$dic$dic,  2))
message("  WAIC: ", round(model_rq1$waic$waic, 2))
message("  Effective parameters (pD):   ", round(model_rq1$dic$p.eff,  1))
message("  Effective parameters (pWAIC):", round(model_rq1$waic$p.eff, 1))

cat("\n--- Fixed effect (intercept = overall log RR) ---\n")
print(model_rq1$summary.fixed)

cat("\n--- Hyperparameters ---\n")
print(model_rq1$summary.hyperpar)

# CPO: conditional predictive ordinate — model adequacy
# Negative sum of log-CPO (NLSCPO): lower is better
cpo_vals <- model_rq1$cpo$cpo
cpo_fail <- model_rq1$cpo$failure
n_fail   <- sum(cpo_fail > 0, na.rm = TRUE)
nlscpo   <- -sum(log(cpo_vals[cpo_fail == 0]), na.rm = TRUE)
message("\n  CPO failures:   ", n_fail, " (", round(100 * n_fail / length(cpo_fail), 2), "%)")
message("  NLSCPO:         ", round(nlscpo, 2))
if (n_fail > 0.05 * length(cpo_fail))
  message("  WARNING: >5% CPO failures — consider inla.cpo() to recompute.")

# =============================================================================
# 7. EXTRACT SPATIAL RANDOM EFFECTS (municipality-level)
# =============================================================================

message("\nExtracting spatial random effects...")

# BYM2 summary: first N entries of summary.random$muni_idx are the
# combined (u_i + v_i) effects; entries N+1..2N are the u_i (structured) only.
N <- n_distinct(inla_df$muni_idx)

spatial_combined <- model_rq1$summary.random$muni_idx[1:N, ] |>
  as.data.frame() |>
  rename(
    muni_idx   = ID,
    s_mean     = mean,
    s_sd       = sd,
    s_q025     = `0.025quant`,
    s_q50      = `0.5quant`,
    s_q975     = `0.975quant`,
    s_mode     = mode
  ) |>
  mutate(
    muni_resid = muni_order[muni_idx],
    # Posterior RR from spatial effect alone: exp(s_mean)
    rr_spatial = exp(s_mean),
    # Exceedance probability: P(spatial RR > 1) = P(s > 0)
    # Approximated from marginals below
    sig_positive = s_q025 > 0,   # 95% CI entirely above 0
    sig_negative = s_q975 < 0    # 95% CI entirely below 0
  )

message("  Municipalities with significantly elevated spatial RR: ",
        sum(spatial_combined$sig_positive))
message("  Municipalities with significantly reduced spatial RR:  ",
        sum(spatial_combined$sig_negative))

# =============================================================================
# 8. EXTRACT TEMPORAL RANDOM EFFECT (RW1)
# =============================================================================

message("Extracting temporal random effect (RW1)...")

temporal <- model_rq1$summary.random$year_idx |>
  as.data.frame() |>
  rename(
    year_idx = ID,
    g_mean   = mean,
    g_sd     = sd,
    g_q025   = `0.025quant`,
    g_q50    = `0.5quant`,
    g_q975   = `0.975quant`,
    g_mode   = mode
  ) |>
  mutate(
    year     = year_idx + 2000L,
    rr_temporal       = exp(g_mean),
    rr_temporal_q025  = exp(g_q025),
    rr_temporal_q975  = exp(g_q975)
  )

cat("\n--- Temporal RW1 posterior (RR scale) ---\n")
print(temporal |> select(year, rr_temporal, rr_temporal_q025, rr_temporal_q975),
      digits = 3)

# =============================================================================
# 9. EXTRACT MUNICIPALITY-YEAR FITTED VALUES (posterior RR)
# =============================================================================

message("Extracting municipality-year fitted relative risks...")

# summary.fitted.values gives posterior mean and quantiles of theta_it = RR_it
fitted_rr <- model_rq1$summary.fitted.values |>
  as.data.frame() |>
  rename(
    rr_mean  = mean,
    rr_sd    = sd,
    rr_q025  = `0.025quant`,
    rr_q50   = `0.5quant`,
    rr_q975  = `0.975quant`
  ) |>
  mutate(
    row_id     = seq_len(n()),
    muni_resid = inla_df$muni_resid[row_id],
    year       = inla_df$year[row_id],
    n_cases    = inla_df$n_cases[row_id],
    E_it       = inla_df$E_it[row_id],
    period     = case_when(
      year <= 2006 ~ "2001-2006",
      year <= 2012 ~ "2007-2012",
      TRUE         ~ "2013-2018"
    ),
    # Exceedance: RR > 1 <=> log(RR) > 0
    excess_prob = ifelse(rr_q025 > 1, TRUE, FALSE)   # simplified flag
  )

message("  Fitted rows: ", nrow(fitted_rr))

# =============================================================================
# 10. EXCEEDANCE PROBABILITIES  P(theta_it > 1)
# =============================================================================
# For each municipality, compute P(RR > 1) from the marginal posterior.
# Uses inla.pmarginal() on the linear predictor marginals (log scale).
# We compute this for the period-averaged spatial effect.

message("Computing exceedance probabilities P(RR > 1) from spatial marginals...")

exceed_prob <- sapply(
  seq_len(N),
  function(i) {
    marg <- model_rq1$marginals.random$muni_idx[[i]]
    if (is.null(marg)) return(NA_real_)
    1 - inla.pmarginal(0, marg)   # P(u_i > 0) = P(spatial RR > 1)
  }
)

spatial_combined$p_exceed <- exceed_prob

message("  Municipalities with P(RR > 1) > 0.80: ",
        sum(exceed_prob > 0.80, na.rm = TRUE))
message("  Municipalities with P(RR > 1) > 0.95: ",
        sum(exceed_prob > 0.95, na.rm = TRUE))

# =============================================================================
# 11. PERIOD-AVERAGED POSTERIOR RR PER MUNICIPALITY
# =============================================================================

message("Computing period-averaged posterior RR per municipality...")

period_rr <- fitted_rr |>
  group_by(muni_resid, period) |>
  summarise(
    rr_period_mean = mean(rr_mean,  na.rm = TRUE),
    rr_period_q025 = mean(rr_q025, na.rm = TRUE),
    rr_period_q975 = mean(rr_q975, na.rm = TRUE),
    total_cases    = sum(n_cases,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(period = factor(period,
                         levels = c("2001-2006", "2007-2012", "2013-2018")))

# =============================================================================
# 12. FIGURES
# =============================================================================

message("\nGenerating figures...")

theme_vl <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(colour = "grey40"),
      panel.grid.minor = element_blank()
    )
}

# ---- Figure 1: Temporal trend (RW1 posterior) ----

p_temporal <- ggplot(temporal, aes(x = year)) +
  geom_ribbon(aes(ymin = rr_temporal_q025, ymax = rr_temporal_q975),
              fill = "#2c7bb6", alpha = 0.25) +
  geom_line(aes(y = rr_temporal), colour = "#2c7bb6", linewidth = 1) +
  geom_point(aes(y = rr_temporal), colour = "#2c7bb6", size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Posterior temporal trend in VL relative risk — Brazil 2001–2018",
    subtitle = "RW1 random effect (exp scale); shaded band = 95% credible interval",
    x        = "Year",
    y        = "Relative risk (RR)"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "rq1_temporal_trend.png"),
       p_temporal, width = 10, height = 5, dpi = 300, bg = "white")
message("  Saved: rq1_temporal_trend.png")

# ---- Figure 2: Posterior mean RR map (period average, 3 panels) ----

make_rr_map <- function(period_label, shp, state_borders, data, rr_cap) {
  pd <- data |> filter(period == period_label)
  shp_pd <- shp |> left_join(pd, by = c("code6" = "muni_resid")) |>
    mutate(rr_plot = pmin(rr_period_mean, rr_cap))
  ggplot(shp_pd) +
    geom_sf(aes(fill = rr_plot), colour = NA, linewidth = 0) +
    geom_sf(data = state_borders, fill = NA,
            colour = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colours  = c("#2c7bb6", "#ffffbf", "#d7191c"),
      na.value = "grey85",
      limits   = c(0, rr_cap),
      values   = scales::rescale(c(0, 1, rr_cap), to = c(0, 1)),
      name     = "Posterior\nmean RR",
      labels   = label_number(accuracy = 0.1),
      guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                 barwidth = 6, barheight = 0.4)
    ) +
    labs(title = period_label) +
    theme_void(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
}

# State borders for overlay
sf_use_s2(FALSE)
state_borders_18 <- muni_2018 |>
  st_make_valid() |>
  mutate(state_code = substr(code6, 1, 2)) |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")
sf_use_s2(TRUE)

rr_cap <- quantile(period_rr$rr_period_mean, 0.99, na.rm = TRUE)

library(patchwork)

p_rr1 <- make_rr_map("2001-2006", muni_2018, state_borders_18, period_rr, rr_cap)
p_rr2 <- make_rr_map("2007-2012", muni_2018, state_borders_18, period_rr, rr_cap)
p_rr3 <- make_rr_map("2013-2018", muni_2018, state_borders_18, period_rr, rr_cap)

map_rr <- (p_rr1 | p_rr2 | p_rr3) +
  plot_annotation(
    title    = "Posterior Mean Relative Risk of VL by Municipality — Brazil 2001–2018",
    subtitle = "INLA Poisson BYM2 + RW1 model; period-averaged posterior RR",
    caption  = "Source: SINAN/IBGE. Agusto Lab, University of Kansas.",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
      plot.caption  = element_text(size = 7,  hjust = 1,   colour = "grey50")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "rq1_rr_map_period.png"),
       map_rr, width = 14, height = 7, dpi = 300, bg = "white")
message("  Saved: rq1_rr_map_period.png")

# ---- Figure 3: Exceedance probability map P(spatial RR > 1) ----

exceed_sf <- muni_2018 |>
  left_join(spatial_combined |> select(muni_resid, p_exceed, rr_spatial),
            by = c("code6" = "muni_resid"))

p_exceed_map <- ggplot(exceed_sf) +
  geom_sf(aes(fill = p_exceed), colour = NA, linewidth = 0) +
  geom_sf(data = state_borders_18, fill = NA,
          colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours  = c("#2c7bb6", "#ffffbf", "#d7191c"),
    na.value = "grey85",
    limits   = c(0, 1),
    name     = "P(spatial RR > 1)",
    labels   = label_percent(accuracy = 1),
    guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 8, barheight = 0.4)
  ) +
  labs(
    title    = "Exceedance Probability: P(Spatial RR > 1) by Municipality",
    subtitle = "Posterior probability of elevated spatial relative risk (BYM2 structured effect)",
    caption  = "Source: SINAN/IBGE. Agusto Lab, University of Kansas."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    plot.caption    = element_text(size = 7, hjust = 1, colour = "grey50"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "rq1_exceedance_map.png"),
       p_exceed_map, width = 8, height = 9, dpi = 300, bg = "white")
message("  Saved: rq1_exceedance_map.png")

# =============================================================================
# 13. SAVE RESULTS
# =============================================================================

message("\nSaving results...")

saveRDS(model_rq1,       file.path(proc_dir, "inla_rq1_model.rds"))
saveRDS(fitted_rr,       file.path(proc_dir, "inla_rq1_results.rds"))
saveRDS(temporal,        file.path(proc_dir, "inla_rq1_temporal.rds"))
saveRDS(spatial_combined, file.path(proc_dir, "inla_rq1_spatial.rds"))
saveRDS(period_rr,       file.path(proc_dir, "inla_rq1_period_rr.rds"))

message("\nSaved to ", proc_dir, ":")
message("  inla_rq1_model.rds      — full INLA model object")
message("  inla_rq1_results.rds    — municipality-year posterior RR")
message("  inla_rq1_temporal.rds   — RW1 temporal trend")
message("  inla_rq1_spatial.rds    — BYM2 spatial effects + exceedance probs")
message("  inla_rq1_period_rr.rds  — period-averaged posterior RR per municipality")

# =============================================================================
# 14. KEY RESULTS SUMMARY FOR THESIS
# =============================================================================

message("\n=== Key results for thesis ===")

# Intercept: overall RR
alpha_mean <- model_rq1$summary.fixed["(Intercept)", "mean"]
message("  Overall intercept (log RR): ", round(alpha_mean, 3),
        "  => RR = ", round(exp(alpha_mean), 3))

# BYM2 hyperparameters
hyp <- model_rq1$summary.hyperpar
message("\n  BYM2 precision (marginal):     ",
        round(hyp["Precision for muni_idx", "mean"], 3))
message("  BYM2 phi (spatial proportion): ",
        round(hyp["Phi for muni_idx", "mean"], 3),
        "  [0=unstructured, 1=fully spatial]")
message("  RW1 precision:                 ",
        round(hyp["Precision for year_idx", "mean"], 3))

# Temporal trend: min and max RR years
message("\n  Peak temporal RR year:   ", temporal$year[which.max(temporal$rr_temporal)],
        "  (RR = ", round(max(temporal$rr_temporal), 3), ")")
message("  Trough temporal RR year: ", temporal$year[which.min(temporal$rr_temporal)],
        "  (RR = ", round(min(temporal$rr_temporal), 3), ")")

# Top municipalities by period-averaged posterior RR
message("\n  Top 10 municipalities by posterior RR (2013–2018 period):")
top10 <- period_rr |>
  filter(period == "2013-2018") |>
  slice_max(rr_period_mean, n = 10) |>
  select(muni_resid, rr_period_mean, rr_period_q025, rr_period_q975, total_cases)
print(top10)

message("\nNext: run 09_inla_rq2.R")

# =============================================================================
# 15. MUNICIPALITY NAMES FOR TOP-10 TABLE
# =============================================================================

source("muni_names.R")

message("\n  Top 10 municipalities by posterior RR (2013–2018) with names:")
top10_named <- period_rr |>
  filter(period == "2013-2018") |>
  slice_max(rr_period_mean, n = 10) |>
  select(muni_resid, rr_period_mean, rr_period_q025, rr_period_q975,
         total_cases) |>
  left_join(muni_lookup, by = c("muni_resid" = "code6")) |>
  mutate(across(where(is.double), ~ round(.x, 2))) |>
  select(muni_name, uf, rr_period_mean, rr_period_q025, rr_period_q975,
         total_cases)

print(top10_named, n = 10)
saveRDS(top10_named, file.path(proc_dir, "inla_rq1_top10_named.rds"))