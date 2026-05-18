# =============================================================================
# Exploratory Data Analysis — VL Brazil (2001–2018)
# Study: Spatial & Spatio-Temporal Dynamics of VL in Brazil (2001–2018)
# Agusto Lab — Mariel Andrianavalondrahona
#
# Run AFTER scripts 01–03 have been executed and their .rds outputs saved.
# This script reads from processed/ and produces no persistent outputs —
# all results are printed or plotted to screen (or saved to figures/).
#
# Sections:
#   1.  Setup & load
#   2.  Case-level descriptives (sinan_confirmed)
#   3.  Annual trends (national)
#   4.  Age and sex profiles
#   5.  HIV co-infection
#   6.  Outcome / CFR (2007–2018)
#   7.  Geographic coverage (municipality-level)
#   8.  Municipality-year incidence distribution
#   9.  Period-level summaries
#  10.  Moran scatter plot (pre-modelling, using raw rates)
#  11.  State-level summaries for context
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)    # for label_comma(), label_percent()
library(sf)        # only needed for section 10 map preview

proc_dir  <- "D:/AIMS/Research phase/R/data/processed"
fig_dir   <- "D:/AIMS/Research phase/R/figures/eda"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: save a ggplot with consistent settings
save_fig <- function(p, name, w = 8, h = 5) {
  ggsave(
    filename = file.path(fig_dir, paste0(name, ".png")),
    plot     = p,
    width    = w, height = h, dpi = 300
  )
  message("  Saved: ", name, ".png")
}

# Shared theme
theme_vl <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "grey40"),
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# 1. LOAD
# =============================================================================

message("Loading processed data...")

sinan <- readRDS(file.path(proc_dir, "sinan_confirmed_2001_2018.rds"))
agg   <- readRDS(file.path(proc_dir, "agg_municipality_year.rds"))

message("  sinan_confirmed: ", nrow(sinan), " rows")
message("  agg (full grid): ", nrow(agg),   " rows")

# Period labels ordered for plots
PERIOD_LEVELS <- c("2001-2006", "2007-2012", "2013-2018")

# =============================================================================
# 2. CASE-LEVEL DESCRIPTIVES
# =============================================================================

message("\n=== 2. Case-level descriptives ===")

cat("\n--- Overall summary ---\n")
cat("Total confirmed cases: ", nrow(sinan), "\n")
cat("Year range:            ", range(sinan$year), "\n")
cat("Unique municipalities: ", length(unique(sinan$muni_resid)), "\n")

cat("\n--- Missing values (key variables) ---\n")
miss <- sinan |>
  summarise(
    age_years_NA  = sum(is.na(age_years)),
    sex_NA        = sum(is.na(sex)),
    hiv_NA        = sum(is.na(hiv)),
    vl_death_NA   = sum(is.na(vl_death))   # expected: all 2001-2006 rows
  )
print(miss)

cat("\n--- Sex distribution ---\n")
sinan |>
  count(sex, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

cat("\n--- HIV status distribution ---\n")
sinan |>
  count(hiv, name = "n") |>
  mutate(
    label = case_when(hiv == 1 ~ "HIV+", hiv == 0 ~ "HIV-", TRUE ~ "Unknown/NA"),
    pct   = round(100 * n / sum(n), 1)
  ) |>
  arrange(hiv) |>
  print()

cat("\n--- Age group distribution ---\n")
sinan |>
  count(age_group, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(age_group) |>
  print()

# =============================================================================
# 3. ANNUAL NATIONAL TRENDS
# =============================================================================

message("\n=== 3. Annual national trends ===")

nat_annual <- agg |>
  group_by(year) |>
  summarise(
    total_cases  = sum(n_cases),
    total_pop    = sum(POPULACAO, na.rm = TRUE),
    total_deaths = sum(n_vl_death, na.rm = TRUE),
    nat_rate     = total_cases / total_pop * 1e5,
    nat_cfr      = ifelse(
      first(year) >= 2007 & total_cases > 0,
      total_deaths / total_cases,
      NA_real_
    ),
    .groups = "drop"
  )

cat("\n--- Annual cases and national incidence rate ---\n")
print(nat_annual |> select(year, total_cases, nat_rate, total_deaths, nat_cfr), n = 18)

# Plot: Annual cases
p_cases <- ggplot(nat_annual, aes(x = year, y = total_cases)) +
  geom_col(fill = "#2c7bb6", colour = "white", linewidth = 0.3) +
  geom_line(colour = "#d7191c", linewidth = 0.8) +
  geom_point(colour = "#d7191c", size = 2) +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title    = "Annual confirmed VL cases — Brazil 2001–2018",
    x        = "Year",
    y        = "Confirmed cases"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_cases, "01_annual_cases")

# Plot: National incidence rate
p_rate <- ggplot(nat_annual, aes(x = year, y = nat_rate)) +
  geom_line(colour = "#d7191c", linewidth = 1) +
  geom_point(colour = "#d7191c", size = 2.5) +
  scale_x_continuous(breaks = 2001:2018) +
  labs(
    title    = "National incidence rate — Brazil 2001–2018",
    subtitle = "Confirmed VL cases per 100,000 inhabitants (annual IBGE denominator)",
    x        = "Year",
    y        = "Incidence per 100,000"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_rate, "02_national_rate")

# Plot: VL deaths and national CFR (2007–2018)
p_cfr <- nat_annual |>
  filter(year >= 2007) |>
  ggplot(aes(x = year)) +
  geom_col(aes(y = total_deaths), fill = "#756bb1", colour = "white",
           linewidth = 0.3) +
  geom_line(aes(y = nat_cfr * max(nat_annual$total_deaths[nat_annual$year >= 2007],
                                  na.rm = TRUE) / max(nat_annual$nat_cfr, na.rm = TRUE)),
            colour = "#d7191c", linewidth = 0.9) +
  scale_x_continuous(breaks = 2007:2018) +
  scale_y_continuous(
    name     = "VL deaths (bars)",
    sec.axis = sec_axis(
      ~ . * max(nat_annual$nat_cfr, na.rm = TRUE) /
        max(nat_annual$total_deaths[nat_annual$year >= 2007], na.rm = TRUE),
      name   = "National CFR (line)",
      labels = label_percent(accuracy = 0.1)
    )
  ) +
  labs(
    title    = "VL deaths and case fatality rate — Brazil 2007–2018",
    subtitle = "Bars = deaths (left axis);  red line = CFR (right axis)",
    x        = "Year"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_cfr, "03_deaths_cfr")

# =============================================================================
# 4. AGE AND SEX PROFILES
# =============================================================================

message("\n=== 4. Age and sex profiles ===")

AGE_LEVELS <- c("<1", "1-4", "5-9", "10-19", "20-39", "40-59", "60+")

# 4a. Age group × sex pyramid
age_sex <- sinan |>
  filter(!is.na(age_group), !is.na(sex)) |>
  count(age_group, sex, name = "n") |>
  mutate(
    age_group = factor(age_group, levels = AGE_LEVELS),
    n_plot    = ifelse(sex == "M", -n, n)    # males left, females right
  )

p_pyramid <- ggplot(age_sex, aes(x = n_plot, y = age_group, fill = sex)) +
  geom_col() +
  scale_x_continuous(
    labels = function(x) label_comma()(abs(x)),
    name   = "Number of cases"
  ) +
  scale_fill_manual(values = c("F" = "#e08080", "M" = "#6080c0"),
                    labels = c("F" = "Female", "M" = "Male")) +
  labs(
    title    = "Age–sex pyramid of confirmed VL cases — Brazil 2001–2018",
    y        = "Age group",
    fill     = NULL
  ) +
  theme_vl()

save_fig(p_pyramid, "04_age_sex_pyramid", w = 7, h = 6)

# 4b. Age distribution over time (by period)
age_period <- sinan |>
  filter(!is.na(age_group)) |>
  mutate(
    period    = case_when(
      year <= 2006 ~ "2001-2006",
      year <= 2012 ~ "2007-2012",
      TRUE         ~ "2013-2018"
    ),
    age_group = factor(age_group, levels = AGE_LEVELS),
    period    = factor(period, levels = PERIOD_LEVELS)
  ) |>
  count(period, age_group, name = "n") |>
  group_by(period) |>
  mutate(pct = n / sum(n)) |>
  ungroup()

p_age_period <- ggplot(age_period, aes(x = age_group, y = pct, fill = period)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  scale_fill_manual(values = c("#2c7bb6", "#fdae61", "#d7191c")) +
  labs(
    title    = "Age distribution of VL cases by period",
    x        = "Age group",
    y        = "Proportion of period cases",
    fill     = "Period"
  ) +
  theme_vl()

save_fig(p_age_period, "05_age_by_period")

# 4c. Sex ratio over time
sex_year <- sinan |>
  filter(!is.na(sex)) |>
  group_by(year, sex) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(year) |>
  mutate(pct = n / sum(n)) |>
  ungroup()

p_sex_year <- ggplot(sex_year, aes(x = year, y = pct, colour = sex)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = c("F" = "#e08080", "M" = "#6080c0"),
                      labels = c("F" = "Female", "M" = "Male")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  labs(
    title  = "Sex distribution of VL cases over time",
    x      = "Year", y = "Proportion", colour = NULL
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_sex_year, "06_sex_over_time")

# =============================================================================
# 5. HIV CO-INFECTION
# =============================================================================

message("\n=== 5. HIV co-infection ===")

hiv_year <- sinan |>
  filter(!is.na(hiv)) |>
  group_by(year) |>
  summarise(
    n_total  = n(),
    n_hiv    = sum(hiv == 1),
    hiv_pct  = n_hiv / n_total,
    .groups  = "drop"
  )

cat("\n--- HIV co-infection rate by year ---\n")
print(hiv_year, n = 18)

p_hiv <- ggplot(hiv_year, aes(x = year, y = hiv_pct)) +
  geom_line(colour = "#756bb1", linewidth = 1) +
  geom_point(colour = "#756bb1", size = 2.5) +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  labs(
    title    = "HIV co-infection rate among confirmed VL cases — Brazil 2001–2018",
    subtitle = "Denominator: cases with known HIV status only",
    x        = "Year",
    y        = "HIV co-infection rate"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_hiv, "07_hiv_coinfection")

# =============================================================================
# 6. OUTCOME / CFR (2007–2018)
# =============================================================================

message("\n=== 6. Outcome / CFR descriptives ===")

outcome_data <- sinan |> filter(period == "2007-2018")

cat("\n--- Outcome distribution (2007–2018 only) ---\n")
outcome_data |>
  mutate(
    outcome_label = case_when(
      outcome_raw == "1" ~ "Cure",
      outcome_raw == "2" ~ "Default/Abandon",
      outcome_raw == "3" ~ "VL death",
      outcome_raw == "4" ~ "Other-cause death",
      outcome_raw == "5" ~ "Transferred",
      TRUE               ~ paste0("Other/Unknown (", outcome_raw, ")")
    )
  ) |>
  count(outcome_label, name = "n") |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(desc(n)) |>
  print()

cat("\n--- CFR by age group (2007–2018) ---\n")
outcome_data |>
  filter(!is.na(age_group)) |>
  group_by(age_group) |>
  summarise(
    cases  = n(),
    deaths = sum(vl_death, na.rm = TRUE),
    cfr    = deaths / cases,
    .groups = "drop"
  ) |>
  mutate(age_group = factor(age_group, levels = AGE_LEVELS)) |>
  arrange(age_group) |>
  print()

cat("\n--- CFR by sex (2007–2018) ---\n")
outcome_data |>
  filter(!is.na(sex)) |>
  group_by(sex) |>
  summarise(
    cases  = n(),
    deaths = sum(vl_death, na.rm = TRUE),
    cfr    = deaths / cases,
    .groups = "drop"
  ) |>
  print()

cat("\n--- CFR by HIV status (2007–2018) ---\n")
outcome_data |>
  filter(!is.na(hiv)) |>
  group_by(hiv) |>
  summarise(
    label  = ifelse(first(hiv) == 1, "HIV+", "HIV-"),
    cases  = n(),
    deaths = sum(vl_death, na.rm = TRUE),
    cfr    = deaths / cases,
    .groups = "drop"
  ) |>
  print()

# Plot: CFR by age group
cfr_age <- outcome_data |>
  filter(!is.na(age_group)) |>
  group_by(age_group) |>
  summarise(
    cases  = n(),
    deaths = sum(vl_death, na.rm = TRUE),
    cfr    = deaths / cases,
    .groups = "drop"
  ) |>
  mutate(age_group = factor(age_group, levels = AGE_LEVELS))

p_cfr_age <- ggplot(cfr_age, aes(x = age_group, y = cfr)) +
  geom_col(fill = "#d7191c") +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  labs(
    title    = "Case fatality rate by age group — Brazil 2007–2018",
    x        = "Age group",
    y        = "CFR"
  ) +
  theme_vl()

save_fig(p_cfr_age, "08_cfr_by_age")

# =============================================================================
# 7. GEOGRAPHIC COVERAGE
# =============================================================================

message("\n=== 7. Geographic coverage ===")

n_munis_ever <- length(unique(sinan$muni_resid))
n_munis_total <- length(unique(agg$muni_resid))

cat("\n--- Municipalities ever reporting >= 1 confirmed case ---\n")
cat("  Ever affected:    ", n_munis_ever, "of", n_munis_total, "total\n")
cat("  Percent affected: ", round(100 * n_munis_ever / n_munis_total, 1), "%\n")

# Number of municipalities affected per year
muni_per_year <- agg |>
  filter(n_cases > 0) |>
  group_by(year) |>
  summarise(n_munis = n_distinct(muni_resid), .groups = "drop")

cat("\n--- Municipalities with >= 1 case by year ---\n")
print(muni_per_year, n = 18)

p_muni <- ggplot(muni_per_year, aes(x = year, y = n_munis)) +
  geom_line(colour = "#2c7bb6", linewidth = 1) +
  geom_point(colour = "#2c7bb6", size = 2.5) +
  scale_x_continuous(breaks = 2001:2018) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title    = "Number of municipalities reporting >= 1 VL case by year",
    subtitle = "Reflects geographic expansion of the epidemic",
    x        = "Year",
    y        = "Municipalities affected"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_muni, "09_municipalities_affected")

# =============================================================================
# 8. MUNICIPALITY-YEAR INCIDENCE DISTRIBUTION
# =============================================================================

message("\n=== 8. Municipality-year incidence distribution ===")

nonzero <- agg |> filter(n_cases > 0, !is.na(incidence_rate))

cat("\n--- Incidence rate among non-zero municipality-years ---\n")
cat("  n (municipality-years): ", nrow(nonzero), "\n")
print(summary(nonzero$incidence_rate))

cat("\n  Percentiles (90th, 95th, 99th, max):\n")
cat("  90th:", round(quantile(nonzero$incidence_rate, 0.90), 2), "\n")
cat("  95th:", round(quantile(nonzero$incidence_rate, 0.95), 2), "\n")
cat("  99th:", round(quantile(nonzero$incidence_rate, 0.99), 2), "\n")
cat("  Max: ", round(max(nonzero$incidence_rate), 2), "\n")

cat("\n  Proportion with rate > 10 per 100k: ",
    round(100 * mean(nonzero$incidence_rate > 10), 1), "%\n")
cat("  Proportion with rate > 50 per 100k: ",
    round(100 * mean(nonzero$incidence_rate > 50), 1), "%\n")
cat("  Proportion with rate > 100 per 100k: ",
    round(100 * mean(nonzero$incidence_rate > 100), 1), "%\n")

# Histogram of incidence (capped at 99th percentile for readability)
cap <- quantile(nonzero$incidence_rate, 0.99)
p_hist <- ggplot(nonzero |> filter(incidence_rate <= cap),
                 aes(x = incidence_rate)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white", linewidth = 0.2) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title    = "Distribution of municipality-year VL incidence rates",
    subtitle = paste0("Non-zero municipality-years only (n = ", nrow(nonzero),
                      "); capped at 99th percentile (", round(cap, 1), " per 100k)"),
    x        = "Incidence per 100,000",
    y        = "Count"
  ) +
  theme_vl()

save_fig(p_hist, "10_incidence_histogram")

# Log-scale for the full distribution
p_log <- ggplot(nonzero, aes(x = log10(incidence_rate + 1))) +
  geom_histogram(bins = 60, fill = "#d7191c", colour = "white", linewidth = 0.2) +
  labs(
    title    = "Distribution of VL incidence (log10 scale)",
    subtitle = "log10(rate + 1); non-zero municipality-years",
    x        = "log10(incidence per 100,000 + 1)",
    y        = "Count"
  ) +
  theme_vl()

save_fig(p_log, "11_incidence_log_histogram")

# =============================================================================
# 9. PERIOD-LEVEL SUMMARIES
# =============================================================================

message("\n=== 9. Period-level summaries ===")

period_summary <- agg |>
  group_by(period) |>
  summarise(
    total_cases   = sum(n_cases),
    total_pop_avg = mean(POPULACAO, na.rm = TRUE) * 6,  # avg annual * 6 years
    period_rate   = total_cases / sum(POPULACAO, na.rm = TRUE) * 1e5,
    n_muni_any    = sum(n_cases > 0),
    total_deaths  = sum(n_vl_death, na.rm = TRUE),
    n_muni_cfr    = sum(n_cases > 0 & year >= 2007),
    .groups       = "drop"
  ) |>
  mutate(period = factor(period, levels = PERIOD_LEVELS)) |>
  arrange(period)

cat("\n--- Period-level summary ---\n")
print(period_summary)

# Period bar chart: cases
p_period_cases <- ggplot(period_summary,
                         aes(x = period, y = total_cases, fill = period)) +
  geom_col(colour = "white") +
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = c("#2c7bb6", "#fdae61", "#d7191c")) +
  labs(
    title = "Total VL cases by analysis period — Brazil",
    x     = "Period",
    y     = "Total confirmed cases"
  ) +
  theme_vl() +
  theme(legend.position = "none")

save_fig(p_period_cases, "12_period_cases", w = 5, h = 5)

# =============================================================================
# 10. MORAN SCATTER PLOT (raw period rate, pre-smoothing)
#     Requires nb lists from script 04. Produces one Moran scatter per period.
#     Note: this is purely descriptive — the formal Moran test is in script 05.
# =============================================================================

message("\n=== 10. Moran scatter plots (requires nb lists) ===")

# Load neighbourhood lists if available
nb_available <- file.exists(
  file.path(proc_dir, "nb_2007.rds"),
  file.path(proc_dir, "nb_2018.rds")
)

if (all(nb_available)) {
  
  library(spdep)
  
  nb_2007 <- readRDS(file.path(proc_dir, "nb_2007.rds"))
  nb_2018 <- readRDS(file.path(proc_dir, "nb_2018.rds"))
  
  shp_dir07 <- "data/shapefiles/2007"
  shp_dir18 <- "data/shapefiles/2018"
  
  shp_path07 <- list.files(shp_dir07, pattern = "\\.shp$",
                           full.names = TRUE, recursive = TRUE)[1]
  shp_path18 <- list.files(shp_dir18, pattern = "\\.shp$",
                           full.names = TRUE, recursive = TRUE)[1]
  
  WATER_BODIES <- c("4300001", "4300002")
  
  muni_2007 <- st_read(shp_path07, quiet = TRUE) |>
    filter(!as.character(GEOCODIG_M) %in% WATER_BODIES) |>
    mutate(code6 = substr(as.character(GEOCODIG_M), 1, 6))
  
  muni_2018 <- st_read(shp_path18, quiet = TRUE) |>
    filter(!as.character(CD_GEOCMU) %in% WATER_BODIES) |>
    mutate(code6 = substr(as.character(CD_GEOCMU), 1, 6))
  
  periods_def <- list(
    list(label = "2001-2006", years = 2001:2006, shp = muni_2007, nb = nb_2007),
    list(label = "2007-2012", years = 2007:2012, shp = muni_2018, nb = nb_2018),
    list(label = "2013-2018", years = 2013:2018, shp = muni_2018, nb = nb_2018)
  )
  
  for (pd in periods_def) {
    
    period_agg <- agg |>
      filter(year %in% pd$years) |>
      group_by(muni_resid) |>
      summarise(
        total_cases = sum(n_cases,   na.rm = TRUE),
        total_pop   = sum(POPULACAO, na.rm = TRUE),
        .groups     = "drop"
      ) |>
      mutate(rate = ifelse(total_pop > 0,
                           total_cases / total_pop * 1e5, 0))
    
    shp_data <- pd$shp |>
      st_drop_geometry() |>
      left_join(period_agg, by = c("code6" = "muni_resid")) |>
      mutate(rate = coalesce(rate, 0))
    
    lw   <- nb2listw(pd$nb, style = "W", zero.policy = TRUE)
    lag  <- lag.listw(lw, shp_data$rate, zero.policy = TRUE)
    
    gm <- moran.test(shp_data$rate, lw,
                     zero.policy = TRUE, alternative = "greater")
    I_val <- round(gm$estimate["Moran I statistic"], 3)
    
    moran_df <- data.frame(rate = shp_data$rate, lag_rate = lag)
    
    p_moran <- ggplot(moran_df, aes(x = scale(rate), y = scale(lag_rate))) +
      geom_point(alpha = 0.25, size = 0.6, colour = "#2c7bb6") +
      geom_smooth(method = "lm", se = FALSE, colour = "#d7191c",
                  linewidth = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      labs(
        title    = paste0("Moran scatter — ", pd$label),
        subtitle = paste0("Global Moran's I = ", I_val,
                          " (p < 0.001)"),
        x        = "Standardised VL incidence rate",
        y        = "Spatially lagged rate"
      ) +
      theme_vl()
    
    save_fig(p_moran, paste0("13_moran_scatter_", gsub("-", "_", pd$label)))
  }
  
} else {
  message("  nb_2007.rds / nb_2018.rds not found — skipping Moran scatter plots.")
  message("  Run 04_bayesian_smoothing.R first to build neighbourhood lists.")
}

# =============================================================================
# 11. STATE-LEVEL SUMMARIES
# =============================================================================

message("\n=== 11. State-level summaries ===")

# State names lookup (two-digit code → UF abbreviation)
state_names <- c(
  "11"="RO","12"="AC","13"="AM","14"="RR","15"="PA",
  "16"="AP","17"="TO","21"="MA","22"="PI","23"="CE",
  "24"="RN","25"="PB","26"="PE","27"="AL","28"="SE",
  "29"="BA","31"="MG","32"="ES","33"="RJ","35"="SP",
  "41"="PR","42"="SC","43"="RS","50"="MS","51"="MT",
  "52"="GO","53"="DF"
)

state_summary <- agg |>
  mutate(uf = state_names[state_code]) |>
  group_by(uf, state_code, period) |>
  summarise(
    total_cases  = sum(n_cases),
    total_pop    = sum(POPULACAO, na.rm = TRUE),
    state_rate   = total_cases / total_pop * 1e5,
    total_deaths = sum(n_vl_death, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(period = factor(period, levels = PERIOD_LEVELS)) |>
  arrange(period, desc(total_cases))

cat("\n--- Top 10 states by cases in each period ---\n")
for (pd in PERIOD_LEVELS) {
  cat("\nPeriod:", pd, "\n")
  state_summary |>
    filter(period == pd) |>
    slice_max(total_cases, n = 10) |>
    select(uf, total_cases, state_rate, total_deaths) |>
    print()
}

# Stacked bar: top 10 states over periods
top_states <- state_summary |>
  group_by(uf) |>
  summarise(grand_total = sum(total_cases), .groups = "drop") |>
  slice_max(grand_total, n = 12) |>
  pull(uf)

p_states <- state_summary |>
  filter(uf %in% top_states) |>
  ggplot(aes(x = reorder(uf, -total_cases),
             y = total_cases,
             fill = period)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = c("#2c7bb6", "#fdae61", "#d7191c")) +
  labs(
    title = "VL cases by state and period — top 12 states",
    x     = "State (UF)",
    y     = "Confirmed cases",
    fill  = "Period"
  ) +
  theme_vl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_states, "14_state_period_cases", w = 10, h = 6)

# =============================================================================
# END
# =============================================================================

message("\n=== EDA complete ===")
message("Figures saved to: ", fig_dir)
message("\nKey findings to note before proceeding to modelling:")
message("  • Check age-group shift between periods (05_age_by_period.png)")
message("  • Check rising HIV co-infection trend (07_hiv_coinfection.png)")
message("  • Note right-skewed incidence distribution (10_incidence_histogram.png)")
message("  • Note municipalities affected over time (09_municipalities_affected.png)")
message("  • Moran scatter confirms positive autocorrelation (13_moran_scatter_*.png)")