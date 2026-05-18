# =============================================================================
# PRIOR SENSITIVITY ANALYSIS — RQ2 INLA Models (Model A + Model B)
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Purpose:
#   Assess whether key inferential conclusions from the RQ2 INLA models
#   (09_inla_rq2_300dpi.R) are robust to the choice of PC priors.
#
#   MODEL A — Municipality-level Binomial BYM2
#     Key conclusions tested:
#       (1) phi = 0.757 — CFR spatial variance is predominantly structured
#       (2) 42 municipalities with P(CFR > national) > 0.95
#       (3) Top-10 high-CFR municipalities (Governador Valadares, Belo Horizonte...)
#
#   MODEL B — Individual-level logistic BYM2 with fixed effects
#     Key conclusions tested:
#       (1) phi = 0.810 — residual spatial variance structured after adjustment
#       (2) Fixed effect ORs (age, sex, HIV, year) stable across priors
#
#   Three alternative prior configurations (same as RQ1 sensitivity):
#
#   M0  Baseline   BYM2: P(sigma > 1) = 0.01, phi P(< 0.5) = 0.50
#   M1  Diffuse    BYM2: P(sigma > 1) = 0.05, phi P(< 0.5) = 0.50
#   M2  Tight      BYM2: P(sigma > 0.5) = 0.01, phi P(< 0.5) = 0.50
#   M3  phi-biased BYM2: P(sigma > 1) = 0.01, phi P(< 0.5) = 0.75
#
# Outputs saved to processed/sensitivity/:
#   sens_rq2A_hyperpar.rds       — Model A hyperparameters across priors
#   sens_rq2A_phi.rds            — Model A phi marginal posteriors
#   sens_rq2A_exceed.rds         — Model A exceedance probability comparison
#   sens_rq2A_top10.rds          — Model A top-10 CFR stability
#   sens_rq2B_hyperpar.rds       — Model B hyperparameters across priors
#   sens_rq2B_phi.rds            — Model B phi marginal posteriors
#   sens_rq2B_fe.rds             — Model B fixed effects (ORs) across priors
#
# Figures saved to figures/sensitivity/:
#   sens_rq2A_hyperpar.png       — Model A hyperparameter CrI comparison
#   sens_rq2A_phi_density.png    — Model A phi marginal posteriors
#   sens_rq2A_exceed_scatter.png — P(CFR > national) baseline vs alternatives
#   sens_rq2A_top10.png          — Top-10 CFR stability dot plot
#   sens_rq2B_hyperpar.png       — Model B hyperparameter CrI comparison
#   sens_rq2B_phi_density.png    — Model B phi marginal posteriors
#   sens_rq2B_fe.png             — Fixed effect ORs across priors (forest plot)
#
# =============================================================================

library(INLA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(scales)
library(patchwork)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
sens_dir <- file.path(proc_dir, "sensitivity")
fig_dir  <- "D:/AIMS/Research phase/R/figures/sensitivity"
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA AND BASELINE MODELS
# =============================================================================

message("Loading data and baseline models...")

sinan <- readRDS(file.path(proc_dir, "sinan_confirmed_2001_2018.rds"))
agg   <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))

shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]
muni_2018  <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

muni_order <- muni_2018$code6
graph_path <- file.path(proc_dir, "nb_2018.graph")
stopifnot("nb_2018.graph not found" = file.exists(graph_path))

# --- Reproduce Model A data (identical to script 09) ---
muni_cfr_data <- agg |>
  filter(year >= 2007) |>
  group_by(muni_resid) |>
  summarise(
    D_i = sum(n_vl_death, na.rm = TRUE),
    C_i = sum(n_cases,    na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(C_i > 0) |>
  mutate(
    crude_cfr = D_i / C_i,
    muni_idx  = match(muni_resid, muni_order)
  ) |>
  filter(!is.na(muni_idx))

nat_cfr       <- sum(muni_cfr_data$D_i) / sum(muni_cfr_data$C_i)
nat_cfr_logit <- log(nat_cfr / (1 - nat_cfr))

# --- Reproduce Model B data (identical to script 09) ---
AGE_REF <- "5-9"

indiv <- sinan |>
  filter(period == "2007-2018", !is.na(vl_death)) |>
  mutate(
    age_group = factor(age_group,
                       levels = c(AGE_REF, "<1", "1-4", "10-19", "20-39", "40-59", "60+")),
    male      = as.integer(sex == "M"),
    hiv_known = as.integer(!is.na(hiv)),
    hiv_pos   = ifelse(is.na(hiv), 0L, hiv),
    year_c    = year - 2012L,
    muni_idx  = match(muni_resid, muni_order),
    death     = vl_death
  ) |>
  filter(!is.na(muni_idx), !is.na(age_group), !is.na(male), !is.na(death))

# Load baseline models
model_A_baseline <- readRDS(file.path(proc_dir, "inla_rq2_modelA.rds"))
model_B_baseline <- readRDS(file.path(proc_dir, "inla_rq2_modelB.rds"))

message("  Baseline models loaded.")
message("  Model A municipalities: ", nrow(muni_cfr_data))
message("  Model B individuals:    ", nrow(indiv))

# Municipality name lookup
source("muni_names.R")   # provides muni_lookup (code6, muni_name, uf)
muni_name_lookup <- muni_lookup |> rename(muni_resid = code6)

# =============================================================================
# 2. DEFINE PRIOR CONFIGURATIONS
# =============================================================================

prior_configs <- list(
  
  M0_baseline = list(
    label      = "M0: Baseline",
    desc       = "BYM2 P(\u03c3>1)=0.01, phi P(<0.5)=0.50",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    )
  ),
  
  M1_diffuse = list(
    label      = "M1: Diffuse",
    desc       = "BYM2 P(\u03c3>1)=0.05, phi P(<0.5)=0.50",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.05)),
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    )
  ),
  
  M2_tight = list(
    label      = "M2: Tight",
    desc       = "BYM2 P(\u03c3>0.5)=0.01, phi P(<0.5)=0.50",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(0.5, 0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 0.5))
    )
  ),
  
  M3_phi_unstructured = list(
    label      = "M3: \u03d5-biased",
    desc       = "BYM2 P(\u03c3>1)=0.01, phi P(<0.5)=0.75",
    prior_bym2 = list(
      prec = list(prior = "pc.prec", param = c(1,   0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 0.75))
    )
  )
)

model_nms  <- names(prior_configs)
model_lbls <- sapply(prior_configs, `[[`, "label")

model_colours <- c(
  "M0: Baseline"  = "#1b7837",
  "M1: Diffuse"   = "#2c7bb6",
  "M2: Tight"     = "#d73027",
  "M3: \u03d5-biased" = "#8856a7"
)

# =============================================================================
# 3. HELPERS
# =============================================================================

fit_modelA <- function(cfg, muni_cfr_data, graph_path) {
  formula <- D_i ~
    1 +
    f(muni_idx,
      model       = "bym2",
      graph       = graph_path,
      scale.model = TRUE,
      constr      = TRUE,
      hyper       = cfg$prior_bym2)
  inla(
    formula  = formula,
    family   = "binomial",
    data     = muni_cfr_data,
    Ntrials  = muni_cfr_data$C_i,
    control.compute = list(dic  = TRUE, waic = TRUE,
                           return.marginals.predictor = TRUE),
    control.predictor = list(compute = TRUE, link = 1),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    verbose = FALSE
  )
}

fit_modelB <- function(cfg, indiv, graph_path) {
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
      hyper       = cfg$prior_bym2)
  inla(
    formula  = formula,
    family   = "binomial",
    data     = indiv,
    Ntrials  = 1,
    control.compute = list(dic  = TRUE, waic = TRUE,
                           return.marginals.predictor = FALSE),
    control.predictor = list(compute = FALSE),
    control.inla = list(strategy = "adaptive", int.strategy = "eb"),
    verbose = FALSE
  )
}

# =============================================================================
# 4. FIT ALTERNATIVE MODELS
# =============================================================================

message("\n--- Fitting Model A alternatives ---")

models_A <- list(M0_baseline = model_A_baseline,
                 M1_diffuse = NULL, M2_tight = NULL, M3_phi_unstructured = NULL)

for (nm in c("M1_diffuse", "M2_tight", "M3_phi_unstructured")) {
  cfg <- prior_configs[[nm]]
  message("  Fitting Model A — ", cfg$label, " ...")
  t0 <- proc.time()
  models_A[[nm]] <- fit_modelA(cfg, muni_cfr_data, graph_path)
  message("    Done in ", round((proc.time() - t0)["elapsed"], 1),
          " sec.  DIC = ", round(models_A[[nm]]$dic$dic, 1))
}

message("\n--- Fitting Model B alternatives ---")
message("  (Model B is individual-level — each run may take 10–25 min)")

models_B <- list(M0_baseline = model_B_baseline,
                 M1_diffuse = NULL, M2_tight = NULL, M3_phi_unstructured = NULL)

for (nm in c("M1_diffuse", "M2_tight", "M3_phi_unstructured")) {
  cfg <- prior_configs[[nm]]
  message("  Fitting Model B — ", cfg$label, " ...")
  t0 <- proc.time()
  models_B[[nm]] <- fit_modelB(cfg, indiv, graph_path)
  message("    Done in ", round((proc.time() - t0)["elapsed"] / 60, 1),
          " min.  DIC = ", round(models_B[[nm]]$dic$dic, 1))
}

message("\nAll models fitted.")

# =============================================================================
# 5. EXTRACT MODEL A RESULTS
# =============================================================================

message("Extracting Model A comparison results...")

N_A <- nrow(muni_cfr_data)

# ---- 5a. Hyperparameters ----

extract_hyperpar <- function(model, model_id, label) {
  hyp <- model$summary.hyperpar
  tibble(
    model_id    = model_id,
    model_label = label,
    parameter   = rownames(hyp),
    mean        = hyp[, "mean"],
    q025        = hyp[, "0.025quant"],
    q975        = hyp[, "0.975quant"]
  )
}

hyperpar_A <- purrr::map2_dfr(
  models_A, model_nms,
  ~ extract_hyperpar(.x, .y, prior_configs[[.y]]$label)
)

# ---- 5b. Phi marginal posteriors ----

phi_A <- purrr::map2(
  models_A, model_nms,
  function(m, nm) {
    marg <- m$marginals.hyperpar[["Phi for muni_idx"]]
    if (is.null(marg)) return(NULL)
    as.data.frame(marg) |>
      mutate(model_id = nm, model_label = prior_configs[[nm]]$label)
  }
) |> bind_rows()

# ---- 5c. Exceedance probabilities P(CFR > national CFR) ----

message("  Computing exceedance probabilities for Model A across all configurations...")

exceed_A <- purrr::map2_dfr(
  models_A, model_nms,
  function(m, nm) {
    p <- sapply(seq_len(N_A), function(i) {
      marg <- m$marginals.fitted.values[[i]]
      if (is.null(marg)) return(NA_real_)
      1 - inla.pmarginal(nat_cfr, marg)
    })
    tibble(
      muni_resid  = muni_cfr_data$muni_resid,
      model_id    = nm,
      model_label = prior_configs[[nm]]$label,
      p_exceed    = p
    )
  }
)

exceed_A_wide <- exceed_A |>
  select(muni_resid, model_id, p_exceed) |>
  pivot_wider(names_from = model_id, values_from = p_exceed)

# ---- 5d. Top-10 CFR stability ----

baseline_top10_cfr <- readRDS(file.path(proc_dir, "inla_rq2_muni_cfr.rds")) |>
  slice_max(cfr_mean, n = 10) |>
  pull(muni_resid)

top10_cfr_sensitivity <- purrr::map2_dfr(
  models_A, model_nms,
  function(m, nm) {
    m$summary.fitted.values |>
      as.data.frame() |>
      rename(cfr_mean = mean, cfr_q025 = `0.025quant`, cfr_q975 = `0.975quant`) |>
      mutate(
        row_id      = seq_len(n()),
        muni_resid  = muni_cfr_data$muni_resid[row_id],
        model_id    = nm,
        model_label = prior_configs[[nm]]$label
      )
  }
) |>
  filter(muni_resid %in% baseline_top10_cfr) |>
  left_join(muni_name_lookup, by = "muni_resid")

# ---- 5e. Model fit table ----

fit_table_A <- purrr::map2_dfr(
  models_A, model_nms,
  function(m, nm) tibble(
    model_id    = nm,
    model_label = prior_configs[[nm]]$label,
    DIC         = round(m$dic$dic,  1),
    WAIC        = round(m$waic$waic, 1),
    phi_mean    = round(m$summary.hyperpar["Phi for muni_idx", "mean"], 3),
    phi_q025    = round(m$summary.hyperpar["Phi for muni_idx", "0.025quant"], 3),
    phi_q975    = round(m$summary.hyperpar["Phi for muni_idx", "0.975quant"], 3),
    n_exceed95  = sum(
      sapply(seq_len(N_A), function(i) {
        marg <- m$marginals.fitted.values[[i]]
        if (is.null(marg)) return(NA_real_)
        1 - inla.pmarginal(nat_cfr, marg)
      }) > 0.95, na.rm = TRUE)
  )
)

cat("\n=== Model A fit comparison ===\n")
print(fit_table_A)

# =============================================================================
# 6. EXTRACT MODEL B RESULTS
# =============================================================================

message("Extracting Model B comparison results...")

N_B <- n_distinct(indiv$muni_idx)

# ---- 6a. Hyperparameters ----

hyperpar_B <- purrr::map2_dfr(
  models_B, model_nms,
  ~ extract_hyperpar(.x, .y, prior_configs[[.y]]$label)
)

# ---- 6b. Phi marginal posteriors ----

phi_B_marg <- purrr::map2(
  models_B, model_nms,
  function(m, nm) {
    marg <- m$marginals.hyperpar[["Phi for muni_idx"]]
    if (is.null(marg)) return(NULL)
    as.data.frame(marg) |>
      mutate(model_id = nm, model_label = prior_configs[[nm]]$label)
  }
) |> bind_rows()

# ---- 6c. Fixed effects (ORs) across priors ----

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

fe_all <- purrr::map2_dfr(
  models_B, model_nms,
  function(m, nm) {
    m$summary.fixed |>
      as.data.frame() |>
      tibble::rownames_to_column("term") |>
      rename(log_or = mean, lo_q025 = `0.025quant`, lo_q975 = `0.975quant`) |>
      mutate(
        or          = exp(log_or),
        or_q025     = exp(lo_q025),
        or_q975     = exp(lo_q975),
        label       = fe_labels[term],
        model_id    = nm,
        model_label = prior_configs[[nm]]$label
      )
  }
)

# ---- 6d. Model B fit table ----

fit_table_B <- purrr::map2_dfr(
  models_B, model_nms,
  function(m, nm) tibble(
    model_id    = nm,
    model_label = prior_configs[[nm]]$label,
    DIC         = round(m$dic$dic,  1),
    WAIC        = round(m$waic$waic, 1),
    phi_mean    = round(m$summary.hyperpar["Phi for muni_idx", "mean"], 3),
    phi_q025    = round(m$summary.hyperpar["Phi for muni_idx", "0.025quant"], 3),
    phi_q975    = round(m$summary.hyperpar["Phi for muni_idx", "0.975quant"], 3)
  )
)

cat("\n=== Model B fit comparison ===\n")
print(fit_table_B)

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

message("Saving sensitivity results...")

saveRDS(hyperpar_A,           file.path(sens_dir, "sens_rq2A_hyperpar.rds"))
saveRDS(phi_A,                file.path(sens_dir, "sens_rq2A_phi.rds"))
saveRDS(exceed_A,             file.path(sens_dir, "sens_rq2A_exceed.rds"))
saveRDS(top10_cfr_sensitivity,file.path(sens_dir, "sens_rq2A_top10.rds"))
saveRDS(fit_table_A,          file.path(sens_dir, "sens_rq2A_fit_table.rds"))
saveRDS(hyperpar_B,           file.path(sens_dir, "sens_rq2B_hyperpar.rds"))
saveRDS(phi_B_marg,           file.path(sens_dir, "sens_rq2B_phi.rds"))
saveRDS(fe_all,               file.path(sens_dir, "sens_rq2B_fe.rds"))
saveRDS(fit_table_B,          file.path(sens_dir, "sens_rq2B_fit_table.rds"))

# =============================================================================
# 8. FIGURES
# =============================================================================

message("Generating sensitivity figures...")

# ---- Figure A1: Model A hyperparameter CrI ----

p_hypA <- hyperpar_A |>
  filter(parameter %in% c("Precision for muni_idx", "Phi for muni_idx")) |>
  mutate(
    param_label = ifelse(parameter == "Phi for muni_idx",
                         "BYM2 \u03d5", "BYM2 precision"),
    model_label = factor(model_label, levels = model_lbls)
  ) |>
  ggplot(aes(x = model_label, y = mean, colour = model_label)) +
  geom_pointrange(aes(ymin = q025, ymax = q975),
                  size = 0.7, linewidth = 0.8) +
  facet_wrap(~ param_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = model_colours, guide = "none") +
  labs(
    title    = "Model A hyperparameters — prior sensitivity",
    subtitle = "Point = posterior mean; bar = 95% credible interval",
    x = NULL, y = "Posterior estimate"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "sens_rq2A_hyperpar.png"),
       p_hypA, width = 9, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_rq2A_hyperpar.png")

# ---- Figure A2: Model A phi density ----

p_phiA <- phi_A |>
  mutate(model_label = factor(model_label, levels = model_lbls)) |>
  ggplot(aes(x = x, y = y, colour = model_label)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Model A: marginal posterior of \u03d5",
    subtitle = "BYM2 spatial mixing parameter — municipality-level CFR",
    x = "\u03d5", y = "Posterior density"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "sens_rq2A_phi_density.png"),
       p_phiA, width = 8, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_rq2A_phi_density.png")

# ---- Figure A3: Exceedance probability scatterplots (Model A) ----

make_exceed_scatter <- function(data, alt_col, alt_label, colour) {
  ggplot(data, aes(x = M0_baseline, y = .data[[alt_col]])) +
    geom_point(alpha = 0.3, size = 0.7, colour = colour) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(limits = c(0, 1), labels = label_percent(accuracy = 1)) +
    scale_y_continuous(limits = c(0, 1), labels = label_percent(accuracy = 1)) +
    labs(title = paste("M0 vs", alt_label),
         x = "P(CFR > national) — M0",
         y = paste("P(CFR > national) —", alt_label)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank())
}

p_excA <- (
  make_exceed_scatter(exceed_A_wide, "M1_diffuse",
                      "M1: Diffuse", model_colours["M1: Diffuse"]) |
    make_exceed_scatter(exceed_A_wide, "M2_tight",
                        "M2: Tight",   model_colours["M2: Tight"]) |
    make_exceed_scatter(exceed_A_wide, "M3_phi_unstructured",
                        "M3: \u03d5-biased", model_colours["M3: \u03d5-biased"])
) +
  plot_annotation(
    title    = "Model A exceedance probability — baseline vs alternative priors",
    subtitle = paste0("Each point = one municipality (n = ", nrow(muni_cfr_data),
                      "). Diagonal = identical inference."),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

ggsave(file.path(fig_dir, "sens_rq2A_exceed_scatter.png"),
       p_excA, width = 13, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_rq2A_exceed_scatter.png")

# ---- Figure A4: Top-10 CFR stability dot plot ----

p_top10A <- top10_cfr_sensitivity |>
  mutate(
    muni_label  = paste0(muni_name, " (", uf, ")"),
    model_label = factor(model_label, levels = model_lbls)
  ) |>
  ggplot(aes(x = cfr_mean * 100, y = reorder(muni_label, cfr_mean),
             colour = model_label)) +
  geom_pointrange(aes(xmin = cfr_q025 * 100, xmax = cfr_q975 * 100),
                  position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 0.7) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Model A top-10 municipality CFR — prior sensitivity",
    subtitle = "Point = posterior mean CFR (%); bar = 95% credible interval",
    x = "Posterior mean CFR (%)", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "sens_rq2A_top10.png"),
       p_top10A, width = 11, height = 7, dpi = 300, bg = "white")
message("  Saved: sens_rq2A_top10.png")

# ---- Figure B1: Model B phi density ----

p_phiB <- phi_B_marg |>
  mutate(model_label = factor(model_label, levels = model_lbls)) |>
  ggplot(aes(x = x, y = y, colour = model_label)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Model B: marginal posterior of \u03d5",
    subtitle = "BYM2 spatial mixing parameter — residual spatial CFR after patient adjustment",
    x = "\u03d5", y = "Posterior density"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "sens_rq2B_phi_density.png"),
       p_phiB, width = 8, height = 5, dpi = 300, bg = "white")
message("  Saved: sens_rq2B_phi_density.png")

# ---- Figure B2: Model B fixed effects across priors ----

p_fe <- fe_all |>
  filter(term != "(Intercept)") |>
  mutate(
    label       = factor(label, levels = rev(unique(label))),
    model_label = factor(model_label, levels = model_lbls)
  ) |>
  ggplot(aes(x = or, y = label, colour = model_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = or_q025, xmax = or_q975),
                  position = position_dodge(width = 0.6),
                  size = 0.4, linewidth = 0.6) +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  scale_colour_manual(values = model_colours, name = "Model") +
  labs(
    title    = "Model B fixed effect odds ratios — prior sensitivity",
    subtitle = "Points = posterior mean OR; bars = 95% credible interval. Log scale.",
    x = "Odds ratio (log scale)", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "sens_rq2B_fe.png"),
       p_fe, width = 10, height = 7, dpi = 300, bg = "white")
message("  Saved: sens_rq2B_fe.png")

# =============================================================================
# 9. CONSOLE SUMMARY FOR THESIS
# =============================================================================

message("\n=== RQ2 prior sensitivity summary ===\n")

cat("--- Model A fit and phi ---\n")
print(fit_table_A |> select(model_label, DIC, phi_mean, phi_q025, phi_q975, n_exceed95))

cat("\n--- Model A: max absolute change in P(CFR > national) vs baseline ---\n")
max_diff_A <- exceed_A_wide |>
  mutate(
    diff_M1 = abs(M1_diffuse          - M0_baseline),
    diff_M2 = abs(M2_tight            - M0_baseline),
    diff_M3 = abs(M3_phi_unstructured - M0_baseline)
  ) |>
  summarise(across(starts_with("diff"), ~ round(max(.x, na.rm = TRUE), 3)))
print(max_diff_A)

cat("\n--- Model B fit and phi ---\n")
print(fit_table_B |> select(model_label, DIC, phi_mean, phi_q025, phi_q975))

cat("\n--- Model B: key fixed effect ORs across priors ---\n")
fe_summary <- fe_all |>
  filter(term %in% c("age_group60+", "male", "hiv_pos", "year_c")) |>
  select(model_label, label, or, or_q025, or_q975) |>
  mutate(across(c(or, or_q025, or_q975), ~ round(.x, 3)))
print(fe_summary)

message("\nSaved to ", sens_dir)
message("Figures saved to ", fig_dir)
message("\nDone. Update chapter 4 with confirmed RQ2 sensitivity numbers.")

# =============================================================================
# END
# =============================================================================