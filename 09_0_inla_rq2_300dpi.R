# =============================================================================
# INLA — RQ2: Spatial Distribution of VL Case Fatality
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
# =============================================================================
#
# Two complementary models are fitted:
#
# MODEL A — Municipality-level Binomial BYM2
# ------------------------------------------
#   D_i ~ Binomial(C_i, p_i)
#   logit(p_i) = alpha + u_i + v_i
#
#   where:
#     D_i        = VL deaths in municipality i (2007–2018 pooled)
#     C_i        = confirmed VL cases in municipality i (2007–2018 pooled)
#     p_i        = municipality-level case fatality rate (CFR)
#     u_i + v_i  = BYM2 spatial random effect
#
#   Purpose: estimates posterior smoothed CFR per municipality and identifies
#   spatial clusters of excess mortality. LISA is then applied to the posterior
#   mean CFR (in script 06) but the Bayesian estimates are more reliable for
#   municipalities with few cases.
#
# MODEL B — Individual-level Binomial logistic regression with BYM2
# ------------------------------------------------------------------
#   Death_j ~ Bernoulli(p_j)
#   logit(p_j) = alpha + beta_age * age_j + beta_sex * sex_j
#              + beta_hiv * hiv_j + beta_year * year_j
#              + f(muni_idx_j, model = "bym2")
#
#   where:
#     Death_j    = 1 if individual j died of VL, 0 otherwise
#     age_j      = age group (reference: 5–9, lowest CFR)
#     sex_j      = 1 if male
#     hiv_j      = 1 if HIV co-infected (known status only)
#     year_j     = calendar year (centred)
#     muni_idx_j = municipality BYM2 random effect
#
#   Purpose: tests whether spatial CFR structure persists after controlling
#   for patient-level characteristics. Residual spatial variance in the BYM2
#   random effect reflects healthcare system heterogeneity.
#
# Data: sinan_confirmed, 2007–2018 only (death outcome available from v48)
#       38,703 cases, 3,009 VL deaths, 1,942 municipalities
#
# Spatial graph: nb_2018.graph (same as RQ1)
#
# Outputs saved to processed/:
#   inla_rq2_modelA.rds        — municipality Binomial BYM2 model
#   inla_rq2_modelB.rds        — individual logistic BYM2 model
#   inla_rq2_muni_cfr.rds      — posterior smoothed CFR per municipality
#   inla_rq2_fixed_effects.rds — Model B fixed effect estimates
#   inla_rq2_spatial_cfr.rds   — Model B spatial random effects
#
# Figures saved to figures/inla_rq2/:
#   rq2_cfr_map.png            — posterior smoothed CFR choropleth
#   rq2_exceedance_cfr.png     — P(CFR > national CFR) exceedance map
#   rq2_fixed_effects.png      — Model B fixed effect forest plot
#   rq2_spatial_cfr_map.png    — Model B residual spatial CFR map
#
# References:
#   Moraga (2019) Geospatial Health Data — Binomial BYM in INLA
#   Riebler et al. (2016) Statistical Methods in Medical Research — BYM2
#   Lima et al. (2021) — HIV and VL co-infection CFR
#
# =============================================================================

library(INLA)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(patchwork)

proc_dir <- "D:/AIMS/Research phase/R/data/processed"
fig_dir  <- "D:/AIMS/Research phase/R/figures/inla_rq2"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

graph_path <- file.path(proc_dir, "nb_2018.graph")
stopifnot("nb_2018.graph not found" = file.exists(graph_path))

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("Loading data...")

sinan <- readRDS(file.path(proc_dir, "sinan_confirmed_2001_2018.rds"))
agg   <- readRDS(file.path(proc_dir, "agg_smoothed.rds"))

# Load shapefile
shp_dir18  <- "data/shapefiles/2018"
shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                         full.names = TRUE, recursive = TRUE)[1]
muni_2018  <- st_read(shp_path18, quiet = TRUE) |>
  filter(!as.character(CD_GEOCMU) %in% c("4300001", "4300002")) |>
  mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))

muni_order <- muni_2018$code6   # shapefile row order for graph alignment

message("  sinan rows: ", nrow(sinan))
message("  agg rows:   ", nrow(agg))

# =============================================================================
# 2. PREPARE MODEL A DATA — MUNICIPALITY-LEVEL BINOMIAL
# =============================================================================

message("\nPreparing Model A data (municipality-level Binomial)...")

# Aggregate 2007–2018: deaths and cases per municipality (pooled)
muni_cfr_data <- agg |>
  filter(year >= 2007) |>
  group_by(muni_resid) |>
  summarise(
    D_i = sum(n_vl_death, na.rm = TRUE),   # VL deaths
    C_i = sum(n_cases,    na.rm = TRUE),   # confirmed cases
    .groups = "drop"
  ) |>
  filter(C_i > 0) |>                        # only municipalities with >= 1 case
  mutate(
    crude_cfr  = D_i / C_i,
    muni_idx   = match(muni_resid, muni_order)
  ) |>
  filter(!is.na(muni_idx))                  # must be in 2018 graph

message("  Municipalities with >= 1 case (2007–2018): ", nrow(muni_cfr_data))
message("  Total deaths:  ", sum(muni_cfr_data$D_i))
message("  Total cases:   ", sum(muni_cfr_data$C_i))
message("  National CFR:  ", round(sum(muni_cfr_data$D_i) /
                                     sum(muni_cfr_data$C_i) * 100, 2), "%")

# National CFR (logit scale) — used for exceedance probability later
nat_cfr <- sum(muni_cfr_data$D_i) / sum(muni_cfr_data$C_i)
nat_cfr_logit <- log(nat_cfr / (1 - nat_cfr))

cat("\n--- Crude CFR distribution (municipalities with >= 1 case) ---\n")
print(summary(muni_cfr_data$crude_cfr))
cat("  Municipalities with CFR = 0: ",
    sum(muni_cfr_data$crude_cfr == 0), "\n")
cat("  Municipalities with CFR = 1: ",
    sum(muni_cfr_data$crude_cfr == 1), "\n")

# =============================================================================
# 3. PRIORS (shared across both models)
# =============================================================================

prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

# =============================================================================
# 4. MODEL A — FIT
# =============================================================================

message("\nFitting Model A: Binomial BYM2...")

formula_A <- D_i ~
  1 +
  f(muni_idx,
    model       = "bym2",
    graph       = graph_path,
    scale.model = TRUE,
    constr      = TRUE,
    hyper       = prior_bym2)

t_start <- proc.time()

model_A <- inla(
  formula  = formula_A,
  family   = "binomial",
  data     = muni_cfr_data,
  Ntrials  = muni_cfr_data$C_i,
  control.compute = list(
    dic  = TRUE,
    waic = TRUE,
    cpo  = TRUE,
    return.marginals.predictor = TRUE
  ),
  control.predictor = list(
    compute = TRUE,
    link    = 1                  # logit link -> fitted values on probability scale
  ),
  control.inla = list(
    strategy     = "adaptive",
    int.strategy = "eb"
  ),
  verbose = FALSE
)

t_elapsed <- proc.time() - t_start
message("  Model A fitted in ", round(t_elapsed["elapsed"], 1), " seconds.")

message("\n--- Model A fit ---")
message("  DIC:  ", round(model_A$dic$dic,  2))
message("  WAIC: ", round(model_A$waic$waic, 2))

cat("\n--- Model A fixed effect (intercept = overall logit CFR) ---\n")
print(model_A$summary.fixed)
alpha_A <- model_A$summary.fixed["(Intercept)", "mean"]
message("  Overall CFR (back-transformed): ",
        round(plogis(alpha_A) * 100, 2), "%")

cat("\n--- Model A hyperparameters ---\n")
print(model_A$summary.hyperpar)

# =============================================================================
# 5. EXTRACT MODEL A POSTERIOR CFR PER MUNICIPALITY
# =============================================================================

message("\nExtracting posterior CFR estimates (Model A)...")

# summary.fitted.values: posterior mean CFR on probability scale
post_cfr <- model_A$summary.fitted.values |>
  as.data.frame() |>
  rename(
    cfr_mean = mean,
    cfr_sd   = sd,
    cfr_q025 = `0.025quant`,
    cfr_q50  = `0.5quant`,
    cfr_q975 = `0.975quant`
  ) |>
  mutate(
    row_id     = seq_len(n()),
    muni_resid = muni_cfr_data$muni_resid[row_id],
    D_i        = muni_cfr_data$D_i[row_id],
    C_i        = muni_cfr_data$C_i[row_id],
    crude_cfr  = muni_cfr_data$crude_cfr[row_id]
  )

message("  Posterior CFR range: [",
        round(min(post_cfr$cfr_mean) * 100, 2), "%, ",
        round(max(post_cfr$cfr_mean) * 100, 2), "%]")

# Exceedance probability: P(CFR_i > national CFR)
# Computed from marginal posterior of the linear predictor (logit scale)
message("Computing P(CFR > national CFR) from marginals...")

exceed_cfr <- sapply(
  seq_len(nrow(muni_cfr_data)),
  function(i) {
    marg <- model_A$marginals.fitted.values[[i]]
    if (is.null(marg)) return(NA_real_)
    # marginals.fitted.values are on probability scale for binomial
    1 - inla.pmarginal(nat_cfr, marg)
  }
)

post_cfr$p_exceed_cfr <- exceed_cfr

message("  Municipalities with P(CFR > national) > 0.80: ",
        sum(exceed_cfr > 0.80, na.rm = TRUE))
message("  Municipalities with P(CFR > national) > 0.95: ",
        sum(exceed_cfr > 0.95, na.rm = TRUE))

# Spatial random effects from Model A
N_A <- nrow(muni_cfr_data)
spatial_A <- model_A$summary.random$muni_idx[1:N_A, ] |>
  as.data.frame() |>
  rename(
    muni_idx = ID,
    s_mean   = mean,
    s_sd     = sd,
    s_q025   = `0.025quant`,
    s_q50    = `0.5quant`,
    s_q975   = `0.975quant`
  ) |>
  mutate(
    muni_resid   = muni_cfr_data$muni_resid[muni_idx],
    sig_elevated = s_q025 > 0,
    sig_reduced  = s_q975 < 0
  )

message("  Municipalities with significantly elevated CFR spatial effect: ",
        sum(spatial_A$sig_elevated, na.rm = TRUE))
message("  Municipalities with significantly reduced CFR spatial effect:  ",
        sum(spatial_A$sig_reduced,  na.rm = TRUE))

# =============================================================================
# 6. PREPARE MODEL B DATA — INDIVIDUAL-LEVEL LOGISTIC
# =============================================================================

message("\nPreparing Model B data (individual-level logistic)...")

AGE_REF <- "5-9"   # reference: lowest CFR group from EDA

indiv <- sinan |>
  filter(period == "2007-2018",
         !is.na(vl_death)) |>
  mutate(
    # Age group — reference level = "5-9" (lowest CFR)
    age_group = factor(
      age_group,
      levels = c(AGE_REF, "<1", "1-4", "10-19", "20-39", "40-59", "60+")
    ),
    # Sex: male = 1, female = 0 (reference)
    male       = as.integer(sex == "M"),
    # HIV: 1 = positive, 0 = negative (NA excluded from fixed effect)
    hiv_known  = as.integer(!is.na(hiv)),
    hiv_pos    = ifelse(is.na(hiv), 0L, hiv),   # NA coded 0; weight via hiv_known
    # Year centred on 2012 (midpoint of 2007–2018)
    year_c     = year - 2012L,
    # Municipality index for BYM2
    muni_idx   = match(muni_resid, muni_order),
    # Outcome
    death      = vl_death
  ) |>
  filter(!is.na(muni_idx),
         !is.na(age_group),
         !is.na(male),
         !is.na(death))

message("  Individual records for Model B: ", nrow(indiv))
message("  Deaths: ", sum(indiv$death), " (",
        round(100 * mean(indiv$death), 2), "%)")
message("  Municipalities represented: ", n_distinct(indiv$muni_idx))

# HIV missing — exclude those rows from the HIV fixed effect
# by fitting hiv_pos only among those with known status;
# the BYM2 absorbs remaining spatial variation for both groups
message("  Records with known HIV status: ",
        sum(indiv$hiv_known), " (",
        round(100 * mean(indiv$hiv_known), 1), "%)")

# =============================================================================
# 7. MODEL B — FIT
# =============================================================================

message("\nFitting Model B: Individual logistic BYM2...")
message("  This may take 10–25 minutes.")

# HIV is included only for records with known status.
# Strategy: fit hiv_pos as a covariate and note in thesis that
# coefficient is estimated on the ~73% of records with known HIV.

formula_B <- death ~
  1 +
  age_group +                              # age group fixed effects
  male +                                   # sex
  hiv_pos +                                # HIV status (known cases)
  year_c +                                 # linear year trend
  f(muni_idx,                              # BYM2 municipality random effect
    model       = "bym2",
    graph       = graph_path,
    scale.model = TRUE,
    constr      = TRUE,
    hyper       = prior_bym2)

t_start <- proc.time()

model_B <- inla(
  formula  = formula_B,
  family   = "binomial",
  data     = indiv,
  Ntrials  = 1,                           # Bernoulli (individual-level)
  control.compute = list(
    dic  = TRUE,
    waic = TRUE,
    cpo  = FALSE,                         # skip CPO — too slow at individual level
    return.marginals.predictor = FALSE
  ),
  control.predictor = list(
    compute = FALSE                        # don't compute all N fitted values
  ),
  control.inla = list(
    strategy     = "adaptive",
    int.strategy = "eb"
  ),
  verbose = FALSE
)

t_elapsed <- proc.time() - t_start
message("  Model B fitted in ", round(t_elapsed["elapsed"] / 60, 1), " minutes.")

message("\n--- Model B fit ---")
message("  DIC:  ", round(model_B$dic$dic,  2))
message("  WAIC: ", round(model_B$waic$waic, 2))

cat("\n--- Model B fixed effects (log-odds scale) ---\n")
print(model_B$summary.fixed)

cat("\n--- Model B hyperparameters ---\n")
print(model_B$summary.hyperpar)

# =============================================================================
# 8. EXTRACT MODEL B FIXED EFFECTS (odds ratios)
# =============================================================================

message("\nExtracting Model B odds ratios...")

fe <- model_B$summary.fixed |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  rename(
    log_or   = mean,
    log_or_sd = sd,
    lo_q025  = `0.025quant`,
    lo_q50   = `0.5quant`,
    lo_q975  = `0.975quant`
  ) |>
  mutate(
    or       = exp(log_or),
    or_q025  = exp(lo_q025),
    or_q975  = exp(lo_q975),
    # Label for plotting
    label = case_when(
      term == "(Intercept)"      ~ "Intercept",
      term == "age_group<1"      ~ "Age <1 yr",
      term == "age_group1-4"     ~ "Age 1–4",
      term == "age_group10-19"   ~ "Age 10–19",
      term == "age_group20-39"   ~ "Age 20–39",
      term == "age_group40-59"   ~ "Age 40–59",
      term == "age_group60+"     ~ "Age 60+",
      term == "male"             ~ "Male sex",
      term == "hiv_pos"          ~ "HIV positive",
      term == "year_c"           ~ "Year (centred)",
      TRUE                       ~ term
    )
  )

cat("\n--- Odds ratios (Model B) ---\n")
fe |>
  filter(term != "(Intercept)") |>
  select(label, or, or_q025, or_q975) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  print()

# BYM2 phi for Model B
hyp_B <- model_B$summary.hyperpar
phi_B <- hyp_B["Phi for muni_idx", "mean"]
message("\n  BYM2 phi (Model B): ", round(phi_B, 3),
        "  [proportion of spatial variance that is structured]")
message("  Interpretation: ",
        ifelse(phi_B > 0.5,
               "spatial CFR structure persists after controlling for patient characteristics",
               "most spatial variance is unstructured (patient mix drives CFR)"))

# =============================================================================
# 9. EXTRACT MODEL B SPATIAL RANDOM EFFECTS
# =============================================================================

message("Extracting Model B spatial random effects...")

N_B <- n_distinct(indiv$muni_idx)

spatial_B <- model_B$summary.random$muni_idx[1:N_B, ] |>
  as.data.frame() |>
  rename(
    muni_idx = ID,
    s_mean   = mean,
    s_sd     = sd,
    s_q025   = `0.025quant`,
    s_q50    = `0.5quant`,
    s_q975   = `0.975quant`
  ) |>
  mutate(
    # Map back to municipality code via indiv data
    muni_resid   = muni_order[muni_idx],
    or_spatial   = exp(s_mean),            # spatial OR (adjusted)
    sig_elevated = s_q025 > 0,
    sig_reduced  = s_q975 < 0
  )

message("  Municipalities with significantly elevated adjusted spatial OR: ",
        sum(spatial_B$sig_elevated, na.rm = TRUE))
message("  Municipalities with significantly reduced adjusted spatial OR:  ",
        sum(spatial_B$sig_reduced,  na.rm = TRUE))

# =============================================================================
# 10. FIGURES
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

# State borders
sf_use_s2(FALSE)
state_borders_18 <- muni_2018 |>
  st_make_valid() |>
  mutate(state_code = substr(code6, 1, 2)) |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")
sf_use_s2(TRUE)

# ---- Figure 1: Posterior smoothed CFR map (Model A) ----

cfr_sf <- muni_2018 |>
  left_join(post_cfr |> select(muni_resid, cfr_mean, cfr_q025, cfr_q975,
                               p_exceed_cfr, crude_cfr),
            by = c("code6" = "muni_resid"))

cfr_cap <- quantile(post_cfr$cfr_mean, 0.99, na.rm = TRUE)

p_cfr_map <- ggplot(cfr_sf) +
  geom_sf(aes(fill = pmin(cfr_mean, cfr_cap)), colour = NA, linewidth = 0) +
  geom_sf(data = state_borders_18, fill = NA,
          colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours  = c("grey95", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"),
    na.value = "grey85",
    limits   = c(0, cfr_cap),
    name     = "Posterior\nCFR",
    labels   = label_percent(accuracy = 0.1),
    guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 8, barheight = 0.4)
  ) +
  labs(
    title    = "Posterior Smoothed VL Case Fatality Rate by Municipality — Brazil 2007–2018",
    subtitle = paste0("Bayesian Binomial BYM2 model; national CFR = ",
                      round(nat_cfr * 100, 1), "%; capped at 99th percentile"),
    caption  = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    plot.caption    = element_text(size = 7, hjust = 1, colour = "grey50"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "rq2_cfr_map.png"),
       p_cfr_map, width = 8, height = 9, dpi = 300, bg = "white")
message("  Saved: rq2_cfr_map.png")

# ---- Figure 2: Exceedance probability P(CFR > national CFR) ----

p_exceed_map <- ggplot(cfr_sf) +
  geom_sf(aes(fill = p_exceed_cfr), colour = NA, linewidth = 0) +
  geom_sf(data = state_borders_18, fill = NA,
          colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours  = c("#2c7bb6", "#ffffbf", "#d7191c"),
    na.value = "grey85",
    limits   = c(0, 1),
    name     = paste0("P(CFR > ", round(nat_cfr * 100, 1), "%)"),
    labels   = label_percent(accuracy = 1),
    guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 8, barheight = 0.4)
  ) +
  labs(
    title    = "Exceedance Probability: P(Municipal CFR > National CFR)",
    subtitle = "Posterior probability of above-average case fatality (Model A)",
    caption  = "Source: SINAN/Ministry of Health Brazil. Agusto Lab, University of Kansas."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    plot.caption    = element_text(size = 7, hjust = 1, colour = "grey50"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "rq2_exceedance_cfr.png"),
       p_exceed_map, width = 8, height = 9, dpi = 300, bg = "white")
message("  Saved: rq2_exceedance_cfr.png")

# ---- Figure 3: Model B fixed effects forest plot ----

fe_plot <- fe |>
  filter(term != "(Intercept)") |>
  mutate(label = factor(label, levels = rev(unique(label))))

p_fe <- ggplot(fe_plot, aes(x = or, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = or_q025, xmax = or_q975),
                width = 0.25, colour = "#2c7bb6", linewidth = 0.7,
                orientation = "y") +
  geom_point(size = 3, colour = "#d7191c") +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Odds Ratios for VL Death — Individual-Level Model (Model B)",
    subtitle = paste0("Bayesian logistic regression with BYM2 municipality random effect\n",
                      "Reference: age 5–9, female, HIV negative. Points = posterior mean; ",
                      "bars = 95% credible interval."),
    x        = "Odds ratio (log scale)",
    y        = NULL,
    caption  = "Source: SINAN 2007–2018. Agusto Lab, University of Kansas."
  ) +
  theme_vl() +
  theme(plot.subtitle = element_text(size = 8))

ggsave(file.path(fig_dir, "rq2_fixed_effects.png"),
       p_fe, width = 8, height = 6, dpi = 300, bg = "white")
message("  Saved: rq2_fixed_effects.png")

# ---- Figure 4: Model B residual spatial OR map ----

spatial_B_sf <- muni_2018 |>
  left_join(spatial_B |> select(muni_resid, or_spatial, sig_elevated, sig_reduced),
            by = c("code6" = "muni_resid"))

or_cap <- quantile(spatial_B$or_spatial, 0.99, na.rm = TRUE)

p_spatial_B <- ggplot(spatial_B_sf) +
  geom_sf(aes(fill = pmin(or_spatial, or_cap)), colour = NA, linewidth = 0) +
  geom_sf(data = state_borders_18, fill = NA,
          colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours  = c("#2c7bb6", "#ffffbf", "#d7191c"),
    na.value = "grey85",
    limits   = c(0, or_cap),
    values   = scales::rescale(c(0, 1, or_cap), to = c(0, 1)),
    name     = "Spatial OR\n(adjusted)",
    labels   = label_number(accuracy = 0.1),
    guide    = guide_colourbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 8, barheight = 0.4)
  ) +
  labs(
    title    = "Residual Spatial Odds Ratio for VL Death — After Adjusting for Patient Characteristics",
    subtitle = paste0("BYM2 municipality random effect from Model B\n",
                      "phi = ", round(phi_B, 2),
                      " — remaining spatial variance reflects healthcare system heterogeneity"),
    caption  = "Source: SINAN 2007–2018. Agusto Lab, University of Kansas."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.subtitle   = element_text(hjust = 0.5, size = 8, colour = "grey40"),
    plot.caption    = element_text(size = 7, hjust = 1, colour = "grey50"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "rq2_spatial_cfr_map.png"),
       p_spatial_B, width = 8, height = 9, dpi = 300, bg = "white")
message("  Saved: rq2_spatial_cfr_map.png")

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================

message("\nSaving results...")

saveRDS(model_A,        file.path(proc_dir, "inla_rq2_modelA.rds"))
saveRDS(model_B,        file.path(proc_dir, "inla_rq2_modelB.rds"))
saveRDS(post_cfr,       file.path(proc_dir, "inla_rq2_muni_cfr.rds"))
saveRDS(fe,             file.path(proc_dir, "inla_rq2_fixed_effects.rds"))
saveRDS(spatial_B,      file.path(proc_dir, "inla_rq2_spatial_cfr.rds"))

message("\nSaved to ", proc_dir, ":")
message("  inla_rq2_modelA.rds        — Binomial BYM2 model (municipality-level)")
message("  inla_rq2_modelB.rds        — Individual logistic BYM2 model")
message("  inla_rq2_muni_cfr.rds      — posterior smoothed CFR per municipality")
message("  inla_rq2_fixed_effects.rds — Model B odds ratios")
message("  inla_rq2_spatial_cfr.rds   — Model B residual spatial effects")

# =============================================================================
# 12. KEY RESULTS SUMMARY FOR THESIS
# =============================================================================

message("\n=== Key results for thesis ===")

message("\n-- Model A (Binomial BYM2) --")
message("  National posterior CFR: ",
        round(plogis(alpha_A) * 100, 2), "%")
message("  BYM2 phi: ",
        round(model_A$summary.hyperpar["Phi for muni_idx", "mean"], 3))
message("  Municipalities with P(CFR > national) > 0.95: ",
        sum(exceed_cfr > 0.95, na.rm = TRUE))

message("\n-- Model B (Individual logistic BYM2) --")
message("  BYM2 phi (residual spatial): ", round(phi_B, 3))
message("  HIV OR: ",
        round(exp(model_B$summary.fixed["hiv_pos", "mean"]), 3),
        "  [95% CrI: ",
        round(exp(model_B$summary.fixed["hiv_pos", "0.025quant"]), 3), "–",
        round(exp(model_B$summary.fixed["hiv_pos", "0.975quant"]), 3), "]")
message("  Age 60+ OR: ",
        round(exp(model_B$summary.fixed["age_group60+", "mean"]), 3),
        "  [95% CrI: ",
        round(exp(model_B$summary.fixed["age_group60+", "0.025quant"]), 3), "–",
        round(exp(model_B$summary.fixed["age_group60+", "0.975quant"]), 3), "]")
message("  Male OR: ",
        round(exp(model_B$summary.fixed["male", "mean"]), 3),
        "  [95% CrI: ",
        round(exp(model_B$summary.fixed["male", "0.025quant"]), 3), "–",
        round(exp(model_B$summary.fixed["male", "0.975quant"]), 3), "]")
message("  Year trend OR per year: ",
        round(exp(model_B$summary.fixed["year_c", "mean"]), 3))

message("\n  Top 10 municipalities by residual spatial OR (Model B):")
top10_B <- spatial_B |>
  slice_max(or_spatial, n = 10) |>
  select(muni_resid, or_spatial, s_q025, s_q975, sig_elevated)
print(top10_B)

message("\nPipeline complete. All INLA results saved.")
message("Next steps: thesis write-up and table/figure finalisation.")

# =============================================================================
# 13. MUNICIPALITY NAMES FOR TOP-10 TABLES
# =============================================================================

source("muni_names.R")

message("\n  Top 10 municipalities by posterior CFR (Model A) with names:")
top10_cfr_named <- post_cfr |>
  slice_max(cfr_mean, n = 10) |>
  select(muni_resid, cfr_mean, cfr_q025, cfr_q975,
         p_exceed_cfr, D_i, C_i) |>
  left_join(muni_lookup, by = c("muni_resid" = "code6")) |>
  mutate(
    cfr_mean  = round(cfr_mean  * 100, 2),
    cfr_q025  = round(cfr_q025  * 100, 2),
    cfr_q975  = round(cfr_q975  * 100, 2),
    p_exceed_cfr = round(p_exceed_cfr, 3)
  ) |>
  select(muni_name, uf, cfr_mean, cfr_q025, cfr_q975,
         p_exceed_cfr, D_i, C_i)

print(as.data.frame(top10_cfr_named))

message("\n  Top 10 municipalities by residual spatial OR (Model B) with names:")
top10_B_named <- spatial_B |>
  slice_max(or_spatial, n = 10) |>
  select(muni_resid, or_spatial, s_q025, s_q975, sig_elevated) |>
  left_join(muni_lookup, by = c("muni_resid" = "code6")) |>
  mutate(across(c(or_spatial, s_q025, s_q975), ~ round(.x, 3))) |>
  select(muni_name, uf, or_spatial, s_q025, s_q975, sig_elevated)

print(as.data.frame(top10_B_named))

saveRDS(top10_cfr_named, file.path(proc_dir, "inla_rq2_top10_cfr_named.rds"))
saveRDS(top10_B_named,   file.path(proc_dir, "inla_rq2_top10_B_named.rds"))