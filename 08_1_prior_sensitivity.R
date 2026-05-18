# =============================================================================
# PRIOR SENSITIVITY ANALYSIS — RQ1 INLA Poisson BYM2 + RW1 Model
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Purpose:
#   Assess whether key inferential conclusions from the baseline INLA model
#   (08_inla_rq1_300dpi.R) are robust to the choice of PC priors.
#
#   Three alternative prior configurations are fitted alongside the baseline:
#
#   M0  Baseline   BYM2: P(sigma > 1)   = 0.01  | RW1: P(sigma > 0.5) = 0.01
#   M1  Diffuse    BYM2: P(sigma > 1)   = 0.05  | RW1: P(sigma > 0.5) = 0.05
#                  (allows more spatial variance and temporal volatility)
#   M2  Tight      BYM2: P(sigma > 0.5) = 0.01  | RW1: P(sigma > 0.3) = 0.01
#                  (stronger shrinkage toward smoother, less clustered surface)
#   M3  phi-prior  BYM2 phi: P(phi < 0.5) = 0.75 (biased toward unstructured)
#                  Other priors same as baseline
#                  (tests sensitivity of phi = 0.764 conclusion)
#
# Outputs saved to processed/sensitivity/:
#   sensitivity_hyperpar.rds     — hyperparameter posteriors for all models
#   sensitivity_phi.rds          — phi posterior comparison
#   sensitivity_temporal.rds     — temporal trend comparison
#   sensitivity_spatial.rds      — spatial effect comparison (mean, p_exceed)
#   sensitivity_top10.rds        — top-10 municipality stability table
#
# Figures saved to figures/sensitivity/:
#   sens_hyperpar.png            — posterior means + CrI for hyperparameters
#   sens_phi_density.png         — marginal posterior of phi across models
#   sens_temporal.png            — RW1 trend overlay (4 models)
#   sens_exceedance_scatter.png  — P(RR > 1): baseline vs alternatives
#   sens_top10.png               — top-10 RR stability dot plot
#
# =============================================================================

library(INLA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(scales)
library(patchwork)

proc_dir  <- "D:/AIMS/Research phase/R/data/processed"
sens_dir  <- file.path(proc_dir, "sensitivity")
fig_dir   <- "D:/AIMS/Research phase/R/figures/sensitivity"
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA AND BASELINE MODEL
# =============================================================================

message("Loading data and baseline model...")

# Re-use the same inla_df built in 08_inla_rq1_300dpi.R.
# To avoid duplicating data prep, we load the saved model and results,
# but we also need inla_df and graph_path to refit alternatives.
# --> Either source the data-prep section of script 08 or reload from saved objects.

agg <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))

shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]
muni_2018  <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

graph_path <- file.path(proc_dir, "nb_2018.graph")
stopifnot("nb_2018.graph not found" = file.exists(graph_path))

national_rate <- sum(agg$n_cases, na.rm = TRUE) /
  sum(agg$POPULACAO, na.rm = TRUE)

muni_order <- muni_2018$code6

inla_df <- agg |>
  mutate(
    E_it      = POPULACAO * national_rate,
    n_cases   = ifelse(is.na(POPULACAO) | POPULACAO <= 0, NA_integer_, n_cases),
    E_it      = ifelse(is.na(E_it) | E_it <= 0, 1e-6, E_it),
    muni_idx  = match(muni_resid, muni_order),
    year_idx  = year - 2000L,
    muni_idx2 = muni_idx
  ) |>
  filter(!is.na(muni_idx))

# Load the baseline model object (already fitted in script 08)
model_baseline <- readRDS(file.path(proc_dir, "inla_rq1_model.rds"))
message("  Baseline model loaded.")

# =============================================================================
# 2. DEFINE PRIOR CONFIGURATIONS
# =============================================================================

# Each entry is a named list with:
#   label      — short model label for plots
#   desc       — one-line description for the summary table
#   prior_bym2 — hyper list passed to f(..., hyper = ...)
#   prior_rw1  — hyper list passed to f(..., hyper = ...)

prior_configs <- list(
  
  M0_baseline = list(
    label      = "M0: Baseline",
    desc       = "BYM2 P(\u03c3>1)=0.01, phi P(<0.5)=0.5, RW1 P(\u03c3>0.5)=0.01",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    ),
    prior_rw1  = list(
      prec = list(prior = "pc.prec", param = c(0.5, 0.01))
    )
  ),
  
  M1_diffuse = list(
    label      = "M1: Diffuse",
    desc       = "BYM2 P(\u03c3>1)=0.05, phi P(<0.5)=0.5, RW1 P(\u03c3>0.5)=0.05",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.05)),  # more permissive
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    ),
    prior_rw1  = list(
      prec = list(prior = "pc.prec", param = c(0.5, 0.05))   # allows more volatility
    )
  ),
  
  M2_tight = list(
    label      = "M2: Tight",
    desc       = "BYM2 P(\u03c3>0.5)=0.01, phi P(<0.5)=0.5, RW1 P(\u03c3>0.3)=0.01",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(0.5, 0.01)),  # stronger shrinkage
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    ),
    prior_rw1  = list(
      prec = list(prior = "pc.prec", param = c(0.3, 0.01))   # smoother trend
    )
  ),
  
  M3_phi_unstructured = list(
    label      = "M3: phi prior\u2192unstructured",
    desc       = "BYM2 P(\u03c3>1)=0.01, phi P(<0.5)=0.75, RW1 P(\u03c3>0.5)=0.01",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 0.75))   # biased toward iid
    ),
    prior_rw1  = list(
      prec = list(prior = "pc.prec", param = c(0.5, 0.01))
    )
  )
)

# =============================================================================
# 3. HELPER: FIT ONE MODEL
# =============================================================================

fit_sensitivity_model <- function(cfg, inla_df, graph_path) {
  
  formula <- n_cases ~
    1 +
    f(muni_idx,
      model       = "bym2",
      graph       = graph_path,
      scale.model = TRUE,
      constr      = TRUE,
      hyper       = cfg$prior_bym2) +
    f(year_idx,
      model       = "rw1",
      scale.model = TRUE,
      constr      = TRUE,
      hyper       = cfg$prior_rw1)
  
  inla(
    formula  = formula,
    family   = "poisson",
    data     = inla_df,
    E        = E_it,
    control.compute = list(
      dic    = TRUE,
      waic   = TRUE,
      config = TRUE,
      return.marginals.predictor = FALSE
    ),
    control.predictor = list(compute = TRUE, link = 1),
    control.inla = list(
      strategy     = "adaptive",
      int.strategy = "eb"
    ),
    verbose = FALSE
  )
}

# =============================================================================
# 4. FIT ALTERNATIVE MODELS  (M0 is already fitted — reuse baseline)
# =============================================================================

message("\nFitting alternative models...")
message("  M0 (baseline): already fitted, reusing.")

models <- list(
  M0_baseline         = model_baseline,
  M1_diffuse          = NULL,
  M2_tight            = NULL,
  M3_phi_unstructured = NULL
)

for (nm in c("M1_diffuse", "M2_tight", "M3_phi_unstructured")) {
  cfg <- prior_configs[[nm]]
  message("  Fitting ", cfg$label, " ...")
  t0 <- proc.time()
  models[[nm]] <- fit_sensitivity_model(cfg, inla_df, graph_path)
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message("    Done in ", elapsed, " min.  DIC = ",
          round(models[[nm]]$dic$dic, 1))
}

message("\nAll models fitted.")

# =============================================================================
# 5. EXTRACT RESULTS FOR COMPARISON
# =============================================================================

message("Extracting comparison results...")

N          <- n_distinct(inla_df$muni_idx)
model_nms  <- names(prior_configs)
model_lbls <- sapply(prior_configs, `[[`, "label")

# ---- 5a. Hyperparameter summary table ----

extract_hyperpar <- function(model, model_id, label) {
  hyp <- model$summary.hyperpar
  tibble(
    model_id   = model_id,
    model_label = label,
    parameter  = rownames(hyp),
    mean       = hyp[, "mean"],
    q025       = hyp[, "0.025quant"],
    q975       = hyp[, "0.975quant"]
  )
}

hyperpar_all <- purrr::map2_dfr(
  models, model_nms,
  ~ extract_hyperpar(.x, .y, prior_configs[[.y]]$label)
)

# Also pull the fixed effect (intercept)
intercepts <- purrr::map2_dfr(
  models, model_nms,
  function(m, nm) {
    fx <- m$summary.fixed
    tibble(
      model_id    = nm,
      model_label = prior_configs[[nm]]$label,
      parameter   = "(Intercept)",
      mean        = fx["(Intercept)", "mean"],
      q025        = fx["(Intercept)", "0.025quant"],
      q975        = fx["(Intercept)", "0.975quant"]
    )
  }
)

hyperpar_all <- bind_rows(hyperpar_all, intercepts)

# Model fit table
fit_table <- purrr::map2_dfr(
  models, model_nms,
  function(m, nm) {
    tibble(
      model_id    = nm,
      model_label = prior_configs[[nm]]$label,
      description = prior_configs[[nm]]$desc,
      DIC         = round(m$dic$dic,  1),
      WAIC        = round(m$waic$waic, 1),
      pD          = round(m$dic$p.eff, 1)
    )
  }
)

cat("\n=== Model fit comparison ===\n")
print(fit_table)

# ---- 5b. Phi marginal posteriors ----

phi_marginals <- purrr::map2(
  models, model_nms,
  function(m, nm) {
    # phi is the second hyperparameter for bym2
    marg <- m$marginals.hyperpar[["Phi for muni_idx"]]
    if (is.null(marg)) return(NULL)
    as.data.frame(marg) |>
      mutate(model_id = nm, model_label = prior_configs[[nm]]$label)
  }
) |> bind_rows()

# ---- 5c. Temporal trend (RW1) ----

extract_temporal <- function(model, model_id, label) {
  model$summary.random$year_idx |>
    as.data.frame() |>
    rename(year_idx = ID, g_mean = mean, g_q025 = `0.025quant`, g_q975 = `0.975quant`) |>
    mutate(
      year            = year_idx + 2000L,
      rr_temporal     = exp(g_mean),
      rr_temporal_q025 = exp(g_q025),
      rr_temporal_q975 = exp(g_q975),
      model_id        = model_id,
      model_label     = label
    )
}

temporal_all <- purrr::map2_dfr(
  models, model_nms,
  ~ extract_temporal(.x, .y, prior_configs[[.y]]$label)
)

# ---- 5d. Spatial random effects + exceedance probabilities ----

extract_spatial <- function(model, model_id, label, N, muni_order) {
  
  sp <- model$summary.random$muni_idx[1:N, ] |>
    as.data.frame() |>
    rename(muni_idx = ID, s_mean = mean, s_q025 = `0.025quant`, s_q975 = `0.975quant`) |>
    mutate(
      muni_resid  = muni_order[muni_idx],
      rr_spatial  = exp(s_mean),
      model_id    = model_id,
      model_label = label
    )
  
  # Exceedance probabilities from marginals
  p_exceed <- sapply(seq_len(N), function(i) {
    marg <- model$marginals.random$muni_idx[[i]]
    if (is.null(marg)) return(NA_real_)
    1 - inla.pmarginal(0, marg)
  })
  
  sp$p_exceed <- p_exceed
  sp
}

message("Computing exceedance probabilities for all models (this takes a few minutes)...")

spatial_all <- purrr::map2_dfr(
  models, model_nms,
  ~ extract_spatial(.x, .y, prior_configs[[.y]]$label, N, muni_order)
)

# ---- 5e. Top-10 municipality stability ----

# Derive the baseline top-10 codes from period_rr (already in memory from the
# baseline fitted values). The saved inla_rq1_top10_named.rds dropped muni_resid
# in its final select(), so we cannot pull codes from it directly.
baseline_top10_codes <- readRDS(file.path(proc_dir, "inla_rq1_period_rr.rds")) |>
  filter(period == "2013-2018") |>
  slice_max(rr_period_mean, n = 10) |>
  pull(muni_resid)

# Build municipality name lookup using the same script as 08_inla_rq1_300dpi.R.
# Provides muni_lookup with columns: code6, muni_name, uf
source("muni_names.R")
muni_name_lookup <- muni_lookup |> rename(muni_resid = code6)

# For each model, extract period-averaged posterior RR for 2013-2018

fitted_all <- purrr::map2_dfr(
  models, model_nms,
  function(m, nm) {
    fv <- m$summary.fitted.values |>
      as.data.frame() |>
      rename(rr_mean = mean, rr_q025 = `0.025quant`, rr_q975 = `0.975quant`) |>
      mutate(
        row_id     = seq_len(n()),
        muni_resid = inla_df$muni_resid[row_id],
        year       = inla_df$year[row_id],
        period     = case_when(
          year <= 2006 ~ "2001-2006",
          year <= 2012 ~ "2007-2012",
          TRUE         ~ "2013-2018"
        ),
        model_id    = nm,
        model_label = prior_configs[[nm]]$label
      )
    fv
  }
)

top10_sensitivity <- fitted_all |>
  filter(period == "2013-2018") |>
  group_by(model_id, model_label, muni_resid) |>
  summarise(
    rr_period_mean = mean(rr_mean, na.rm = TRUE),
    rr_period_q025 = mean(rr_q025, na.rm = TRUE),
    rr_period_q975 = mean(rr_q975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(model_id) |>
  mutate(rank = rank(-rr_period_mean)) |>
  ungroup() |>
  filter(muni_resid %in% baseline_top10_codes) |>
  left_join(muni_name_lookup, by = "muni_resid")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

message("Saving sensitivity results...")

saveRDS(hyperpar_all,      file.path(sens_dir, "sensitivity_hyperpar.rds"))
saveRDS(phi_marginals,     file.path(sens_dir, "sensitivity_phi.rds"))
saveRDS(temporal_all,      file.path(sens_dir, "sensitivity_temporal.rds"))
saveRDS(spatial_all,       file.path(sens_dir, "sensitivity_spatial.rds"))
saveRDS(top10_sensitivity, file.path(sens_dir, "sensitivity_top10.rds"))
saveRDS(fit_table,         file.path(sens_dir, "sensitivity_fit_table.rds"))

# =============================================================================
# 7. FIGURES
# =============================================================================

message("Generating sensitivity figures...")

model_colours <- c(
  "M0: Baseline"                    = "#1b7837",
  "M1: Diffuse"                     = "#2c7bb6",
  "M2: Tight"                       = "#d73027",
  "M3: phi prior\u2192unstructured" = "#8856a7"
)

# ---- Figure 1: Hyperparameter posterior means + CrI ----

# Focus on the three key hyperparameters
key_params <- c(
  "Precision for muni_idx" = "BYM2 precision",
  "Phi for muni_idx"       = "BYM2 \u03d5 (spatial proportion)",
  "Precision for year_idx" = "RW1 precision"
)

p_hyp <- hyperpar_all |>
  filter(parameter %in% names(key_params)) |>
  mutate(
    param_label  = key_params[parameter],
    model_label  = factor(model_label, levels = model_lbls)
  ) |>
  ggplot(aes(x = model_label, y = mean, colour = model_label)) +
  geom_pointrange(aes(ymin = q025, ymax = q975),
                  size = 0.7, linewidth = 0.8,
                  position = position_dodge(width = 0.3)) +
  facet_wrap(~ param_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = model_colours, guide = "none") +
  labs(
    title    = "Hyperparameter posterior estimates across prior configurations",
    subtitle = "Point = posterior mean; bar = 95% credible interval",
    x        = NULL,
    y        = "Posterior mean (± 95% CrI)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 9),
    plot.title       = element_text(face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "sens_hyperpar.png"),
       p_hyp, width = 12, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_hyperpar.png")

# ---- Figure 2: Phi marginal posterior densities ----

p_phi <- phi_marginals |>
  mutate(model_label = factor(model_label, levels = model_lbls)) |>
  ggplot(aes(x = x, y = y, colour = model_label)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Marginal posterior of \u03d5 (BYM2 spatial mixing parameter)",
    subtitle = "Higher \u03d5 = more spatially structured clustering",
    x        = "\u03d5",
    y        = "Posterior density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "sens_phi_density.png"),
       p_phi, width = 8, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_phi_density.png")

# ---- Figure 3: Temporal trend overlay ----

p_temporal <- temporal_all |>
  mutate(model_label = factor(model_label, levels = model_lbls)) |>
  ggplot(aes(x = year, y = rr_temporal,
             colour = model_label, fill = model_label)) +
  geom_ribbon(aes(ymin = rr_temporal_q025, ymax = rr_temporal_q975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = model_colours, name = "Model") +
  scale_fill_manual(  values = model_colours, name = "Model") +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Posterior temporal trend (RW1) — prior sensitivity",
    subtitle = "RR on exponential scale; shaded band = 95% credible interval (M0 only for clarity)",
    x        = "Year",
    y        = "Relative risk (RR)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(fig_dir, "sens_temporal.png"),
       p_temporal, width = 11, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_temporal.png")

# ---- Figure 4: Exceedance probability scatter — baseline vs alternatives ----

exceed_wide <- spatial_all |>
  select(muni_resid, model_id, p_exceed) |>
  pivot_wider(names_from = model_id, values_from = p_exceed)

make_exceed_scatter <- function(data, alt_col, alt_label, colour) {
  ggplot(data, aes(x = M0_baseline, y = .data[[alt_col]])) +
    geom_point(alpha = 0.25, size = 0.6, colour = colour) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(limits = c(0, 1), labels = label_percent(accuracy = 1)) +
    scale_y_continuous(limits = c(0, 1), labels = label_percent(accuracy = 1)) +
    labs(
      title = paste("M0 vs", alt_label),
      x     = "P(RR > 1) — M0 Baseline",
      y     = paste("P(RR > 1) —", alt_label)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
}

p_exc1 <- make_exceed_scatter(exceed_wide, "M1_diffuse",
                              "M1: Diffuse",       model_colours["M1: Diffuse"])
p_exc2 <- make_exceed_scatter(exceed_wide, "M2_tight",
                              "M2: Tight",         model_colours["M2: Tight"])
p_exc3 <- make_exceed_scatter(exceed_wide, "M3_phi_unstructured",
                              "M3: phi\u2192unstr", model_colours["M3: phi prior\u2192unstructured"])

p_scatter <- (p_exc1 | p_exc2 | p_exc3) +
  plot_annotation(
    title    = "Exceedance probability P(spatial RR > 1) — baseline vs alternative priors",
    subtitle = "Each point = one municipality (n = 5,570). Points on diagonal = identical inference.",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave(file.path(fig_dir, "sens_exceedance_scatter.png"),
       p_scatter, width = 13, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_exceedance_scatter.png")

# ---- Figure 5: Top-10 municipality RR stability dot plot ----

p_top10 <- top10_sensitivity |>
  mutate(
    muni_label  = paste0(muni_name, " (", uf, ")"),
    model_label = factor(model_label, levels = model_lbls)
  ) |>
  ggplot(aes(x = rr_period_mean, y = reorder(muni_label, rr_period_mean),
             colour = model_label)) +
  geom_pointrange(aes(xmin = rr_period_q025, xmax = rr_period_q975),
                  position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 0.7) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Top-10 municipality posterior RR (2013\u20132018) — prior sensitivity",
    subtitle = "Point = posterior mean; bar = 95% credible interval",
    x        = "Posterior mean relative risk",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(fig_dir, "sens_top10.png"),
       p_top10, width = 11, height = 7, dpi = 300, bg = "white")
message("  Saved: sens_top10.png")

# =============================================================================
# 8. CONSOLE SUMMARY FOR THESIS
# =============================================================================

message("\n=== Prior sensitivity summary ===\n")

cat("--- Model fit ---\n")
print(fit_table)

cat("\n--- Phi (spatial mixing parameter) across models ---\n")
phi_summary <- hyperpar_all |>
  filter(parameter == "Phi for muni_idx") |>
  select(model_label, mean, q025, q975) |>
  rename(`Posterior mean` = mean, `2.5%` = q025, `97.5%` = q975)
print(phi_summary)

cat("\n--- Exceedance threshold counts (P > 0.95) ---\n")
exceed_counts <- spatial_all |>
  group_by(model_label) |>
  summarise(
    `P(RR>1) > 0.95` = sum(p_exceed > 0.95, na.rm = TRUE),
    `P(RR>1) > 0.80` = sum(p_exceed > 0.80, na.rm = TRUE),
    .groups = "drop"
  )
print(exceed_counts)

cat("\n--- Peak and trough temporal RR years ---\n")
temporal_summary <- temporal_all |>
  group_by(model_label) |>
  summarise(
    peak_year  = year[which.max(rr_temporal)],
    peak_RR    = round(max(rr_temporal), 3),
    trough_year = year[which.min(rr_temporal)],
    trough_RR  = round(min(rr_temporal), 3),
    .groups    = "drop"
  )
print(temporal_summary)

cat("\n--- Max absolute change in P(RR>1) vs baseline (by municipality) ---\n")
max_diff <- exceed_wide |>
  mutate(
    diff_M1 = abs(M1_diffuse          - M0_baseline),
    diff_M2 = abs(M2_tight            - M0_baseline),
    diff_M3 = abs(M3_phi_unstructured - M0_baseline)
  ) |>
  summarise(across(starts_with("diff"), ~ round(max(.x, na.rm = TRUE), 3)))
print(max_diff)

message("\nSaved to ", sens_dir)
message("Figures saved to ", fig_dir)
message("\nDone. Proceed to report robustness in Section 3.1 of the research essay.")

# =============================================================================
# END
# =============================================================================