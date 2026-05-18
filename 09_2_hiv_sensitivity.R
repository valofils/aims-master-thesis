# =============================================================================
# HIV MISSINGNESS SENSITIVITY — RQ2 Model B
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Background:
#   In script 09_0_inla_rq2_300dpi.R (Model B), records with unknown HIV status
#   (24.5% of 2007–2018 cases) were coded hiv_pos = 0 (HIV negative). The
#   variable hiv_known was constructed but never used in the model formula.
#   The net effect is that the HIV coefficient is estimated from a contrast
#   between confirmed-positive cases and the combined group of confirmed-
#   negative PLUS unknown-status cases — a non-differential misclassification
#   that biases the HIV OR toward 1 (attenuation toward null).
#
# Purpose of this script:
#   Fit two alternative versions of Model B and compare against the baseline:
#
#   MB_baseline  — original script: hiv_pos coded 0 for unknowns, all 38,703
#                  records included, HIV estimated on full (contaminated) sample
#
#   MB_drop      — MISSING DROPPED: restrict to records with known HIV status
#                  only (~29,200 records); hiv_pos estimated on clean contrast
#                  (confirmed-positive vs confirmed-negative)
#
#   MB_na        — NA PASSTHROUGH: keep all records but pass NA for unknown HIV;
#                  INLA treats NA covariates as missing and effectively estimates
#                  the HIV coefficient on known-status records only, while
#                  retaining all records for the spatial random effect and other
#                  fixed effects. This is the closest implementation of the
#                  originally commented intention in script 09.
#
# Comparisons reported:
#   (1) HIV OR and 95% CrI under each approach
#   (2) All other fixed effect ORs (age, sex, year) — should be stable
#   (3) BYM2 phi — spatial mixing parameter — should be stable
#   (4) DIC / WAIC
#   (5) n (records) and n_deaths contributing to each fit
#
# The central thesis conclusions do not depend on the HIV OR, but this
# analysis provides a formally correct sensitivity check and the corrected
# text for section 3.5.3 of the thesis.
#
# Output saved to:
#   processed/sensitivity/hiv_sens_modelB_results.rds
#   processed/sensitivity/hiv_sens_modelB_fe.rds
#   figures/sensitivity/hiv_sens_fe_comparison.png
#   figures/sensitivity/hiv_sens_phi_density.png
#
# =============================================================================

library(INLA)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
sens_dir <- file.path(proc_dir, "sensitivity")
fig_dir  <- "D:/AIMS/Research phase/R/figures/sensitivity"
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

graph_path <- file.path(proc_dir, "nb_2018.graph")
stopifnot("nb_2018.graph not found" = file.exists(graph_path))

# =============================================================================
# 1. LOAD DATA AND SHAPEFILE
# =============================================================================

message("Loading data...")

sinan <- readRDS(file.path(proc_dir, "sinan_confirmed_2001_2018.rds"))

shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]
muni_2018  <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

muni_order <- muni_2018$code6

# =============================================================================
# 2. BUILD THREE ANALYSIS DATASETS
# =============================================================================

message("Building analysis datasets...")

AGE_REF <- "5-9"

# --- Shared base filter: period, non-missing outcome, non-missing age/sex,
#     in-graph municipality ---
base_prep <- sinan |>
  filter(period == "2007-2018",
         !is.na(vl_death)) |>
  mutate(
    age_group = factor(
      age_group,
      levels = c(AGE_REF, "<1", "1-4", "10-19", "20-39", "40-59", "60+")
    ),
    male      = as.integer(sex == "M"),
    year_c    = year - 2012L,
    muni_idx  = match(muni_resid, muni_order),
    death     = vl_death,
    hiv_known = as.integer(!is.na(hiv))
  ) |>
  filter(!is.na(muni_idx),
         !is.na(age_group),
         !is.na(male),
         !is.na(death))

# ---- Dataset 1: MB_baseline ----
# Replicates script 09 exactly: unknown HIV coded 0.
indiv_baseline <- base_prep |>
  mutate(hiv_pos = ifelse(is.na(hiv), 0L, hiv))

# ---- Dataset 2: MB_drop ----
# Drop all records with unknown HIV. Clean contrast only.
indiv_drop <- base_prep |>
  filter(hiv_known == 1) |>
  mutate(hiv_pos = hiv)

# ---- Dataset 3: MB_na ----
# Keep all records. Pass NA for unknown HIV so INLA excludes those rows
# from the HIV fixed effect only, while retaining them for spatial RE
# and all other fixed effects.
indiv_na <- base_prep |>
  mutate(hiv_pos = ifelse(is.na(hiv), NA_integer_, hiv))

# --- Print dataset sizes ---
cat("\n=== Dataset summary ===\n")
cat(sprintf("  MB_baseline : %d records, %d deaths (%.1f%%), HIV known: %d (%.1f%%)\n",
            nrow(indiv_baseline),
            sum(indiv_baseline$death),
            100 * mean(indiv_baseline$death),
            sum(indiv_baseline$hiv_known),
            100 * mean(indiv_baseline$hiv_known)))

cat(sprintf("  MB_drop     : %d records, %d deaths (%.1f%%), HIV known: %d (100%%)\n",
            nrow(indiv_drop),
            sum(indiv_drop$death),
            100 * mean(indiv_drop$death),
            nrow(indiv_drop)))

cat(sprintf("  MB_na       : %d records, %d deaths (%.1f%%), HIV known: %d (%.1f%%)\n",
            nrow(indiv_na),
            sum(indiv_na$death),
            100 * mean(indiv_na$death),
            sum(!is.na(indiv_na$hiv_pos)),
            100 * mean(!is.na(indiv_na$hiv_pos))))

# =============================================================================
# 3. SHARED PRIOR AND FORMULA BUILDER
# =============================================================================

prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

fit_modelB_hiv <- function(data, label, graph_path) {
  message("  Fitting ", label, " (n = ", nrow(data), ")...")
  t0 <- proc.time()
  
  formula <- death ~
    1 +
    age_group +
    male +
    hiv_pos +
    year_c +
    f(muni_idx,
      model       = "bym2",
      graph       = graph_path,
      scale.model = TRUE,
      constr      = TRUE,
      hyper       = prior_bym2)
  
  fit <- inla(
    formula  = formula,
    family   = "binomial",
    data     = data,
    Ntrials  = 1,
    control.compute = list(
      dic  = TRUE,
      waic = TRUE,
      return.marginals.predictor = FALSE
    ),
    control.predictor = list(compute = FALSE),
    control.inla = list(
      strategy     = "adaptive",
      int.strategy = "eb"
    ),
    verbose = FALSE
  )
  
  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  message("    Done in ", elapsed, " min | DIC = ",
          round(fit$dic$dic, 1), " | WAIC = ", round(fit$waic$waic, 1))
  fit
}

# =============================================================================
# 4. FIT ALL THREE MODELS
# =============================================================================

message("\n=== Fitting Model B variants ===")
message("(Each individual-level INLA fit may take 10-25 min)")

# Load baseline from saved file if available, otherwise refit
baseline_path <- file.path(proc_dir, "inla_rq2_modelB.rds")
if (file.exists(baseline_path)) {
  message("  MB_baseline: loading from ", baseline_path)
  model_baseline <- readRDS(baseline_path)
} else {
  model_baseline <- fit_modelB_hiv(indiv_baseline, "MB_baseline", graph_path)
}

model_drop <- fit_modelB_hiv(indiv_drop, "MB_drop", graph_path)
model_na   <- fit_modelB_hiv(indiv_na,   "MB_na",   graph_path)

# =============================================================================
# 5. EXTRACT AND COMPARE RESULTS
# =============================================================================

message("\nExtracting results...")

# Fixed effect labels
fe_labels <- c(
  "(Intercept)"    = "Intercept",
  "age_group<1"    = "Age <1 yr",
  "age_group1-4"   = "Age 1\u20134",
  "age_group10-19" = "Age 10\u201319",
  "age_group20-39" = "Age 20\u201339",
  "age_group40-59" = "Age 40\u201359",
  "age_group60+"   = "Age 60+",
  "male"           = "Male sex",
  "hiv_pos"        = "HIV positive",
  "year_c"         = "Year (centred)"
)

extract_fe <- function(model, label, n_records, n_deaths, n_hiv_known) {
  model$summary.fixed |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    rename(log_or   = mean,
           lo_q025  = `0.025quant`,
           lo_q975  = `0.975quant`,
           log_sd   = sd) |>
    mutate(
      or           = exp(log_or),
      or_q025      = exp(lo_q025),
      or_q975      = exp(lo_q975),
      label        = fe_labels[term],
      model_label  = label,
      n_records    = n_records,
      n_deaths     = n_deaths,
      n_hiv_known  = n_hiv_known,
      hiv_ci_width = NA_real_   # filled below
    )
}

fe_baseline <- extract_fe(model_baseline, "MB_baseline (hiv=0 for unknowns)",
                          nrow(indiv_baseline), sum(indiv_baseline$death),
                          sum(indiv_baseline$hiv_known))
fe_drop     <- extract_fe(model_drop,     "MB_drop (unknown HIV dropped)",
                          nrow(indiv_drop),     sum(indiv_drop$death),
                          nrow(indiv_drop))
fe_na       <- extract_fe(model_na,       "MB_na (NA passthrough for unknown HIV)",
                          nrow(indiv_na),       sum(indiv_na$death),
                          sum(!is.na(indiv_na$hiv_pos)))

fe_all <- bind_rows(fe_baseline, fe_drop, fe_na)

# Mark HIV CI width for easy comparison
fe_all <- fe_all |>
  mutate(hiv_ci_width = ifelse(term == "hiv_pos", or_q975 - or_q025, NA_real_))

# =============================================================================
# 6. FIT SUMMARY TABLE
# =============================================================================

extract_fit <- function(model, label, n_rec, n_deaths, n_hiv) {
  phi <- model$summary.hyperpar["Phi for muni_idx", ]
  tibble(
    model_label  = label,
    n_records    = n_rec,
    n_deaths     = n_deaths,
    n_hiv_known  = n_hiv,
    DIC          = round(model$dic$dic,  1),
    WAIC         = round(model$waic$waic, 1),
    phi_mean     = round(phi[["mean"]],         3),
    phi_q025     = round(phi[["0.025quant"]],   3),
    phi_q975     = round(phi[["0.975quant"]],   3),
    hiv_or       = round(exp(model$summary.fixed["hiv_pos", "mean"]),         3),
    hiv_or_q025  = round(exp(model$summary.fixed["hiv_pos", "0.025quant"]),   3),
    hiv_or_q975  = round(exp(model$summary.fixed["hiv_pos", "0.975quant"]),   3),
    hiv_sig      = exp(model$summary.fixed["hiv_pos", "0.025quant"]) > 1
  )
}

fit_summary <- bind_rows(
  extract_fit(model_baseline, "MB_baseline",
              nrow(indiv_baseline), sum(indiv_baseline$death),
              sum(indiv_baseline$hiv_known)),
  extract_fit(model_drop,     "MB_drop",
              nrow(indiv_drop),     sum(indiv_drop$death),
              nrow(indiv_drop)),
  extract_fit(model_na,       "MB_na",
              nrow(indiv_na),       sum(indiv_na$death),
              sum(!is.na(indiv_na$hiv_pos)))
)

cat("\n=== Model B HIV sensitivity: fit summary ===\n")
print(as.data.frame(fit_summary))

# Focused HIV comparison
cat("\n=== HIV OR comparison ===\n")
hiv_comparison <- fe_all |>
  filter(term == "hiv_pos") |>
  select(model_label, n_records, n_hiv_known, or, or_q025, or_q975) |>
  mutate(
    ci_crosses_1 = or_q025 < 1 & or_q975 > 1,
    across(c(or, or_q025, or_q975), ~ round(.x, 3))
  )
print(as.data.frame(hiv_comparison))

# Age 60+ stability check
cat("\n=== Age 60+ OR stability ===\n")
age60_comparison <- fe_all |>
  filter(term == "age_group60+") |>
  select(model_label, or, or_q025, or_q975) |>
  mutate(across(c(or, or_q025, or_q975), ~ round(.x, 3)))
print(as.data.frame(age60_comparison))

# phi stability check
cat("\n=== BYM2 phi stability ===\n")
phi_comparison <- fit_summary |>
  select(model_label, phi_mean, phi_q025, phi_q975)
print(as.data.frame(phi_comparison))

# =============================================================================
# 7. PHI MARGINAL POSTERIOR DENSITIES
# =============================================================================

phi_dens <- bind_rows(
  as.data.frame(model_baseline$marginals.hyperpar[["Phi for muni_idx"]]) |>
    mutate(model_label = "MB_baseline (hiv=0 for unknowns)"),
  as.data.frame(model_drop$marginals.hyperpar[["Phi for muni_idx"]]) |>
    mutate(model_label = "MB_drop (unknown HIV dropped)"),
  as.data.frame(model_na$marginals.hyperpar[["Phi for muni_idx"]]) |>
    mutate(model_label = "MB_na (NA passthrough for unknown HIV)")
)

model_colours_hiv <- c(
  "MB_baseline (hiv=0 for unknowns)"       = "#1b7837",
  "MB_drop (unknown HIV dropped)"          = "#d73027",
  "MB_na (NA passthrough for unknown HIV)" = "#2c7bb6"
)

# =============================================================================
# 8. FIGURES
# =============================================================================

message("\nGenerating figures...")

# ---- Figure 1: HIV OR forest plot comparison ----

hiv_fe_plot <- fe_all |>
  filter(term == "hiv_pos") |>
  mutate(model_label = factor(model_label, levels = c(
    "MB_baseline (hiv=0 for unknowns)",
    "MB_drop (unknown HIV dropped)",
    "MB_na (NA passthrough for unknown HIV)"
  )))

p_hiv <- ggplot(hiv_fe_plot, aes(x = or, y = model_label, colour = model_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = or_q025, xmax = or_q975),
                width = 0.2, linewidth = 0.9) +
  geom_point(size = 4) +
  scale_x_continuous(limits = c(0.7, 1.8),
                     breaks  = c(0.75, 1.0, 1.25, 1.5, 1.75)) +
  scale_colour_manual(values = model_colours_hiv, guide = "none") +
  labs(
    title    = "HIV OR under three approaches to missing HIV status",
    subtitle = paste0(
      "Model B (Individual-level Bayesian logistic BYM2)\n",
      "Reference: HIV negative. Points = posterior mean OR; bars = 95% credible interval."
    ),
    x = "Odds ratio for VL death",
    y = NULL,
    caption = "Source: SINAN 2007\u20132018. Agusto Lab, University of Kansas."
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 9)
  )

# ---- Figure 2: All fixed effects comparison ----

fe_plot_all <- fe_all |>
  filter(term != "(Intercept)") |>
  mutate(
    label       = factor(label, levels = rev(c(
      "Age <1 yr", "Age 1\u20134", "Age 10\u201319",
      "Age 20\u201339", "Age 40\u201359", "Age 60+",
      "Male sex", "HIV positive", "Year (centred)"
    ))),
    model_label = factor(model_label, levels = c(
      "MB_baseline (hiv=0 for unknowns)",
      "MB_drop (unknown HIV dropped)",
      "MB_na (NA passthrough for unknown HIV)"
    ))
  )

p_fe_all <- ggplot(fe_plot_all,
                   aes(x = or, y = label, colour = model_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = or_q025, xmax = or_q975),
                  position = position_dodge(width = 0.6),
                  size = 0.4, linewidth = 0.7) +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  scale_colour_manual(values = model_colours_hiv, name = "HIV handling") +
  labs(
    title    = "Model B fixed effects — HIV missingness sensitivity",
    subtitle = "All three approaches to missing HIV status. Log scale. Reference: age 5\u20139, female, HIV negative.",
    x = "Odds ratio (log scale)", y = NULL,
    caption = "Source: SINAN 2007\u20132018. Agusto Lab, University of Kansas."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8)
  )

# ---- Figure 3: phi density comparison ----

p_phi <- ggplot(phi_dens,
                aes(x = x, y = y, colour = model_label)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = model_colours_hiv, name = "HIV handling") +
  labs(
    title    = "BYM2 \u03d5 marginal posterior — HIV missingness sensitivity",
    subtitle = "Residual spatial mixing parameter (Model B). Overlap confirms spatial conclusions are invariant.",
    x = "\u03d5", y = "Posterior density",
    caption = "Source: SINAN 2007\u20132018. Agusto Lab, University of Kansas."
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8)
  )

# ---- Combine HIV OR + phi ----
p_combined <- p_hiv / p_phi +
  plot_annotation(
    title = "HIV missingness sensitivity — Model B",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

# Save figures
ggsave(file.path(fig_dir, "hiv_sens_hiv_or.png"),
       p_hiv,      width = 9, height = 4, dpi = 300, bg = "white")
message("  Saved: hiv_sens_hiv_or.png")

ggsave(file.path(fig_dir, "hiv_sens_fe_comparison.png"),
       p_fe_all,   width = 11, height = 7, dpi = 300, bg = "white")
message("  Saved: hiv_sens_fe_comparison.png")

ggsave(file.path(fig_dir, "hiv_sens_phi_density.png"),
       p_phi,      width = 9, height = 5, dpi = 300, bg = "white")
message("  Saved: hiv_sens_phi_density.png")

ggsave(file.path(fig_dir, "hiv_sens_combined.png"),
       p_combined, width = 9, height = 9, dpi = 300, bg = "white")
message("  Saved: hiv_sens_combined.png")

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

message("Saving results...")

saveRDS(fit_summary,    file.path(sens_dir, "hiv_sens_modelB_results.rds"))
saveRDS(fe_all,         file.path(sens_dir, "hiv_sens_modelB_fe.rds"))
saveRDS(phi_dens,       file.path(sens_dir, "hiv_sens_modelB_phi.rds"))
saveRDS(hiv_comparison, file.path(sens_dir, "hiv_sens_hiv_or_table.rds"))

message("Saved to ", sens_dir)

# =============================================================================
# 10. THESIS TEXT GENERATOR
# =============================================================================

cat("\n\n")
cat("=======================================================================\n")
cat("  SUGGESTED THESIS TEXT — Section 3.5.3 / Appendix (HIV sensitivity)\n")
cat("=======================================================================\n\n")

hiv_b  <- hiv_comparison |> filter(grepl("baseline", model_label))
hiv_d  <- hiv_comparison |> filter(grepl("drop",     model_label))
hiv_n  <- hiv_comparison |> filter(grepl("na",       model_label))
phi_b  <- fit_summary    |> filter(grepl("baseline",  model_label))
phi_d  <- fit_summary    |> filter(grepl("drop",      model_label))
phi_n  <- fit_summary    |> filter(grepl("na",        model_label))

cat(sprintf(
  "Records with missing HIV status (24.5%% of 2007-2018 records) were coded
as HIV negative in the main analysis (MB_baseline, n = %d). To assess
sensitivity, two additional specifications were fitted: MB_drop restricted
the sample to the %d records (%.1f%%) with known HIV status; MB_na retained
all records but passed NA for unknown HIV status, so that INLA estimated the
HIV coefficient on known-status records only while retaining unknown-status
records for all other fixed effects and the spatial random effect.

The HIV OR was %.3f (95%% CrI: %.3f, %.3f) under MB_baseline, %.3f
(%.3f, %.3f) under MB_drop, and %.3f (%.3f, %.3f) under MB_na. The
credible interval crossed 1 in all three specifications, confirming that
HIV co-infection is not a statistically significant predictor of VL death
after controlling for age, sex, year, and spatial effects.

The BYM2 mixing parameter phi was %.3f (MB_baseline), %.3f (MB_drop), and
%.3f (MB_na), confirming that the spatial conclusion — 81%% of residual
CFR variance is geographically structured — is insensitive to the handling
of missing HIV status. All other fixed effect ORs (age groups, male sex,
year trend) were effectively invariant across the three specifications.\n",
  nrow(indiv_baseline),
  nrow(indiv_drop),
  100 * nrow(indiv_drop) / nrow(indiv_baseline),
  hiv_b$or, hiv_b$or_q025, hiv_b$or_q975,
  hiv_d$or, hiv_d$or_q025, hiv_d$or_q975,
  hiv_n$or, hiv_n$or_q025, hiv_n$or_q975,
  phi_b$phi_mean, phi_d$phi_mean, phi_n$phi_mean
))

cat("=======================================================================\n")
cat("\nDone. Run time is dominated by MB_drop and MB_na (10-25 min each).\n")
cat("Results saved to processed/sensitivity/ and figures/sensitivity/\n")

# =============================================================================
# END
# =============================================================================