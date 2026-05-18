# Spatial and Spatio-Temporal Dynamics of Visceral Leishmaniasis in Brazil (2001–2018)

**Author:** Mariel Andrianavalondrahona  
**Supervisor:** Dr. Folashade B. Agusto (Agusto Lab, University of Kansas)  
**Institution:** African Institute for Mathematical Sciences (AIMS), Ghana  
**Programme:** Master's in Mathematical Sciences, Research Essay 2026

---

## Overview

This repository contains the R analysis pipeline for a spatial epidemiology study of visceral leishmaniasis (VL) in Brazil over 2001–2018. The study addresses two research questions:

- **RQ1:** What are the spatiotemporal incidence patterns of VL across Brazilian municipalities (2001–2018)?
- **RQ2:** What are the spatial distribution and individual-level predictors of VL case fatality (2007–2018)?

Both questions are answered using Bayesian spatial models fitted with R-INLA (BYM2 specification).

---

## Data Sources

The following data are **not included** in this repository and must be obtained separately:

| Dataset | Source | Notes |
|---|---|---|
| SINAN VL notification records | Brazilian Ministry of Health (DATASUS) | Two CSV files covering 2001–2006 and 2007–2018 |
| Municipal population estimates | IBGE (POPTBR files) | ZIP archives containing CSV (2001–2013) or DBF (2014–2018) |
| Municipal boundary shapefiles | IBGE | Downloaded automatically by `00_Download_shapefiles.R` |

Place the SINAN CSVs and IBGE population ZIPs in `data/DISTRIB_VI_LEISHMANIASIS/` and `data/Pop_Brazil_2001_2018/` respectively, then update the `data_dir` and `pop_dir` path variables in scripts `01` and `02` to match your local paths.

---

## Repository Structure

```
.
├── data/
│   ├── shapefiles/          # Downloaded by 00_Download_shapefiles.R
│   │   ├── 2007/
│   │   └── 2018/
│   ├── DISTRIB_VI_LEISHMANIASIS/   # SINAN raw data (not included)
│   ├── Pop_Brazil_2001_2018/       # IBGE population files (not included)
│   └── processed/                  # Intermediate .rds files written by scripts
├── figures/                 # All PNG/PDF outputs at 300 DPI
└── outputs/
    └── maps/                # Publication-quality map files
```

---

## Script Pipeline

Run scripts in numbered order. Each script reads `.rds` files saved by the previous step from `data/processed/`.

### Data preparation

| Script | Purpose | Key outputs |
|---|---|---|
| `00_Download_shapefiles.R` | Downloads 2007 and 2018 IBGE municipal boundary shapefiles from IBGE FTP; constructs 6-digit municipality code (`code6`); removes water-body polygons (Lagoa Mirim, Lagoa dos Patos) | `data/shapefiles/2007/`, `data/shapefiles/2018/` |
| `01_load_sinan.R` | Reads and harmonises SINAN VL records for 2001–2006 and 2007–2018; parses heterogeneous age encoding; standardises municipality codes to 6 digits (×10 from SINAN 5-digit); flags confirmed cases, VL deaths, and HIV status | `processed/sinan_confirmed_2001_2018.rds` |
| `02_load_population.R` | Reads 18 annual IBGE population ZIP archives (CSV format 2001–2013, DBF format 2014–2018); normalises all municipality codes to 6 digits | `processed/population_2001_2018.rds` |
| `03_aggregate.R` | Joins cases to population; constructs complete municipality × year panel (including zero-case rows); computes crude incidence rates per 100,000 and CFR (2007–2018 only) | `processed/agg_municipality_year.rds` |
| `muni_names.R` | Builds a `code6 → municipality name + state` lookup table from the 2018 shapefile attribute table | In-memory object `muni_lookup` (sourced by other scripts) |

### Exploratory analysis

| Script | Purpose | Key outputs |
|---|---|---|
| `04_bayesian_smoothing.R` | Applies global empirical Bayes smoothing (EBest, Clayton & Kaldor 1987) per year; constructs queen-contiguity spatial neighbour lists for 2007 and 2018 shapefiles using `poly2nb()`; exports INLA graph files | `processed/agg_smoothed.rds`, `processed/nb_2007.rds`, `processed/nb_2018.rds`, `data/processed/nb_2018.graph` |
| `05_eda.R` | Produces descriptive statistics and figures: annual case counts, incidence trends, age/sex profiles, HIV co-infection rates, CFR trends, geographic spread (number of affected municipalities), Moran scatter plot | `figures/eda/*.png` |
| `06_lisa.R` | Computes Global Moran's I and LISA (Local Moran's I via `localmoran_perm()`) for three six-year periods (2001–06, 2007–12, 2013–18) on smoothed incidence rates | `processed/map_lisa_2001_2006.rds`, `processed/map_lisa_2007_2012.rds`, `processed/map_lisa_2013_2018.rds`, `processed/lisa_annual.rds` |
| `07_maps.R` | Assembles publication-quality choropleth and LISA cluster maps for the three analysis periods; spotlight maps for 2001, 2009, 2018 | `outputs/maps/*.png` |
| `brazil_map.R` | Produces a Brazil locator map with South America inset using `rnaturalearth` | `brazil_locator_map.pdf` |

### INLA modelling — RQ1 (incidence)

| Script | Purpose | Key outputs |
|---|---|---|
| `08_0_inla_rq1_300dpi.R` | Fits a Poisson BYM2 + RW1 spatio-temporal model: `Y_it ~ Poisson(E_it · θ_it)`, `log(θ_it) = α + BYM2_i + RW1_t`; computes posterior relative risk (RR) and exceedance probabilities P(RR > 1) per municipality-year | `processed/inla_rq1_model.rds`, `processed/inla_rq1_results.rds`, `processed/inla_rq1_temporal.rds`, `processed/inla_rq1_spatial.rds`, `figures/inla_rq1/*.png` |
| `08_1_prior_sensitivity.R` | Refits the RQ1 model under four PC-prior configurations (baseline, diffuse, tight, phi-biased) to assess robustness of φ, temporal trend, and top-10 municipality rankings | `processed/sensitivity/sensitivity_*.rds`, `figures/sensitivity/sens_*.png` |

### INLA modelling — RQ2 (case fatality)

| Script | Purpose | Key outputs |
|---|---|---|
| `09_0_inla_rq2_300dpi.R` | Fits two complementary models. **Model A:** municipality-level Binomial BYM2 for pooled 2007–2018 CFR. **Model B:** individual-level Bernoulli BYM2 with fixed effects for age group, sex, HIV status, and calendar year | `processed/inla_rq2_modelA.rds`, `processed/inla_rq2_modelB.rds`, `processed/inla_rq2_muni_cfr.rds`, `processed/inla_rq2_fixed_effects.rds`, `figures/inla_rq2/*.png` |
| `09_1_prior_sensitivity_rq2.R` | Prior sensitivity for both RQ2 models (four configurations); assesses stability of φ, exceedance probabilities, top-10 CFR municipalities, and fixed effect ORs | `processed/sensitivity/sens_rq2*.rds`, `figures/sensitivity/sens_rq2*.png` |
| `09_2_hiv_sensitivity.R` | Sensitivity analysis for HIV missingness in Model B: compares three approaches (recode unknowns as negative; drop unknowns; pass `NA` to INLA) | `processed/sensitivity/hiv_sens_modelB_*.rds`, `figures/sensitivity/hiv_sens_*.png` |

### LISA for CFR

| Script | Purpose | Key outputs |
|---|---|---|
| `10_lisa_cfr.R` | LISA on pooled 2007–2018 CFR (consistent with INLA RQ2 structure); Global Moran's I and cluster classification | `processed/lisa_cfr_pooled.rds`, `processed/map_lisa_cfr_pooled.rds` |
| `10_lisa_cfr_maps.R` | Publication-quality LISA cluster map and Moran scatter plot for CFR | `figures/lisa_cfr/*.png` |

---

## Key Methodological Notes

**Municipality codes.** SINAN uses 5-digit codes; IBGE population files and shapefiles use 6- or 7-digit codes. Conversion: SINAN code × 10 = 6-digit code (e.g. 24081 → 240810). Shapefiles use the first 6 digits of their 7-digit code.

**Water-body exclusion.** Two non-municipality polygon rows (codes 4300001 Lagoa Mirim, 4300002 Lagoa dos Patos) appear in both shapefile editions and are removed before any spatial operations.

**Shapefile geometry fix.** The 2018 shapefile contains one invalid polygon (Pará, state 15). All scripts using this file call `sf_use_s2(FALSE)` and `st_make_valid()` before `poly2nb()` or `st_union()`.

**Global vs local EB smoothing.** Global EB (`EBest`, spdep) is used instead of local EB because approximately 87% of municipalities report zero cases in any given year. When a municipality has cases but all neighbours report zero, local EB over-shrinks estimates to near zero. The global prior (national mean rate) is more appropriate for this sparse, spatially concentrated disease.

**RQ2 restricted to 2007–2018.** VL death outcome codes are absent from SINAN files for 2001–2006; CFR analysis is therefore limited to the 12-year period where the outcome is structurally recorded.

**LISA permutation testing.** spdep ≥ 1.4.2 deprecates the `nsim` argument in `localmoran()`; permutation-based p-values are obtained via `localmoran_perm()`.

**Three six-year periods.** Period boundaries (2001–06, 2007–12, 2013–18) follow the documented epidemiological phases in Bruhn et al. (2024), de Melo et al. (2023), and Lima et al. (2021). Six-year aggregation stabilises rates before LISA.

---

## R Package Requirements

```r
install.packages(c(
  "sf", "dplyr", "tidyr", "readr", "ggplot2", "patchwork",
  "scales", "RColorBrewer", "spdep", "spatialreg", "foreign",
  "cowplot"
))

# R-INLA (not on CRAN)
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)

# rnaturalearth family (for brazil_map.R)
install.packages(c("rnaturalearth", "rnaturalearthdata"))
install.packages(
  "rnaturalearthhires",
  repos = "https://ropensci.r-universe.dev"
)
```

---

## References

- Clayton, D. & Kaldor, J. (1987). Empirical Bayes estimates of age-standardized relative risks for use in disease mapping. *Biometrics*, 43(3), 671–681.
- Riebler, A. et al. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. *Statistical Methods in Medical Research*, 25(4), 1145–1165.
- Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models. *Journal of the Royal Statistical Society B*, 71(2), 319–392.
- Bruhn, F.R.P. et al. (2024). Spatial distribution of visceral leishmaniasis in Brazil.
- de Melo, G.C. et al. (2023). Temporal trends and spatial distribution of visceral leishmaniasis in Brazil.
- Lima, I.D. et al. (2021). Epidemiology of visceral leishmaniasis in Brazil, 2001–2018.
