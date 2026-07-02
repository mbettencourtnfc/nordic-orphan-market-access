# Nordic Orphan Drug Market Access

R-based analysis of time-to-reimbursement and cross-country decision concordance for 51 EMA-approved orphan drugs evaluated by **TLV** (Sweden), **Medicinrådet** (Denmark), and **Nye Metoder** (Norway) between 2018–2026.

---

## Research questions

1. **Time to access** — After EMA marketing authorisation, how long do Nordic HTA bodies take to make a positive reimbursement decision? Does this differ by country, technology type, or evidence base?
2. **Cross-country concordance** — When the same drug is evaluated in all three countries, do they reach the same decision? What predicts disagreement?
3. **Predictors of a positive outcome** — Controlling for drug-level heterogeneity, which country and technology factors are associated with approval?

---

## Dataset

**Coverage:** 51 EMA-approved orphan drugs × 3 countries = 174 decision records (2018–2026).  
**Source:** Manually curated from published HTA body decision documents — TLV reimbursement decisions, Medicinrådet recommendations, and Nye Metoder assessments. All documents are publicly available.

**Key variables:**
- Drug name, INN, therapeutic area, disease
- Technology group (small molecule, monoclonal antibody, gene therapy, RNA-targeted, enzyme replacement therapy, recombinant protein, cell therapy, other)
- Severity tier (ultra-severe → moderate; ~70% missing — coded from decision text)
- Evidence type (RCT, single-arm, real-world/registry/natural history)
- HTA outcome (positive / restricted positive / negative)
- Decision date (available for ~136/174 records; remainder approximated as 1 July of the decision year)
- EMA first marketing authorisation date (merged from EMA medicines data table)
- Sweden ICER ranges in SEK/QALY where reported

**Honest caveats:**
- N=51 drugs is small for statistical modelling — all models should be interpreted as exploratory
- Severity tier is missing for ~70% of records (not systematically coded in HTA documents)
- Decision dates use a July-1 approximation for ~38 records where only the year was available
- Negative decisions are treated as censored in survival analysis (a simplification — resubmission is possible)
- Manual coding of Evidence_Type introduces subjectivity

---

## Methods

### `00_fetch_ema_ma_dates.R`
Downloads the EMA medicines data table, extracts the first marketing authorisation date per drug (not amendment dates), applies fuzzy matching for combination INNs, and saves `data/raw/ema_ma_dates.csv`. Coverage: 49/51 drugs (2 remain NA).

### `01_data_cleaning.R`
Loads the master dataset, fixes column types, merges EMA dates, computes `time_to_decision_years`, flags records with negative or implausibly long times, and saves `data/nordic_orphan_clean.csv`.

### `02_time_to_access.R` — Survival analysis
- **Kaplan-Meier curves** stratified by country: cumulative probability of having received a positive decision by years from EMA authorisation
- **Log-rank test**: overall and pairwise (Bonferroni-corrected) comparison of country-specific curves
- **Cox proportional hazards** (two models): clustered standard errors by drug to account for the same drug appearing in multiple countries
- HR > 1 = faster rate of positive decision; reference: Norway, Cell therapy, RCT evidence

### `03_concordance.R` — Concordance and mixed-effects logistic regression
- Reshapes to wide format (one row per drug), classifies each drug as concordant-all-positive, concordant-all-negative, or discordant
- **Mixed-effects logistic regression** with drug random intercept `(1|Drug)` — absorbs unmeasured drug-level factors (clinical value, price, evidence quality) so fixed effects estimate country/technology contributions net of those
- OR > 1 = higher odds of positive decision vs reference (Norway, Cell therapy)

---

## Key findings

### Time to access
| Country | N | Positive decisions | Median time to approval (yr) |
|---------|---|-------------------|------------------------------|
| Norway  | 45 | 16 (36%) | 6.61 |
| Denmark | 37 | 23 (62%) | 2.28 |
| Sweden  | 47 | 25 (53%) | 5.39 |

- Log-rank p=0.13 overall (underpowered at N=51 drugs). Pairwise: Norway vs Denmark p=0.18.
- Cox Model A (N=129 valid records): Denmark HR=2.85\*\*\*, Sweden HR=1.85\* vs Norway after adjusting for technology type and evidence.
- **RNA-targeted drugs** reach a positive decision fastest (HR=6.98\*\*\* vs cell therapy).
- **Real-world/registry evidence** is associated with faster approval (HR=8.62\*) — counter-intuitive, but reflects that orphan diseases often cannot run RCTs, and HTA bodies accept natural history data for high-severity indications.
- **Cell therapy** has the lowest approval rates: 0% in Norway and Denmark, 43% in Sweden.

### Cross-country concordance
- All 51 drugs were evaluated in all 3 countries.
- **49% of drugs (25/51) generated discordant decisions** across countries.
- Only **22% (11/51)** received unanimous approval; **29% (15/51)** were unanimously rejected.
- Dominant discordance pattern: Norway rejects drugs that Denmark and Sweden approve.
- Notable exception: Sweden rejected Elfabrio and Fabrazyme (enzyme replacement therapies) that both Norway and Denmark approved — Sweden's stricter cost-effectiveness threshold in play.

### Predictors of approval (mixed-effects logistic regression)
Reference: Norway, Cell therapy. N=171 records, 50 drugs, drug random intercept.

| Predictor | OR | p |
|-----------|-----|---|
| Denmark vs Norway | 1.82 | 0.23 |
| Sweden vs Norway | 2.05 | 0.15 |
| Small molecule vs Cell therapy | 24.9 | 0.018* |
| RNA-targeted | 14.3 | 0.067. |
| Enzyme replacement therapy | 10.3 | 0.092. |
| Gene therapy | 2.3 | 0.55 |

- Country effects trend in the expected direction but do not reach significance — the model is underpowered.
- Technology type is the strongest predictor: small molecules have 25× higher odds of approval vs cell therapies.
- Drug random-effect SD=1.51 — most outcome variation is at the drug level (intrinsic evidence quality, price, unmet need), not explained by country or technology type alone.

---

## Repo structure

```
nordic-orphan-market-access/
├── data/
│   ├── raw/
│   │   └── ema_ma_dates.csv          # EMA first MA dates (output of script 00)
│   └── nordic_orphan_clean.csv       # clean analysis dataset
├── notebooks/
│   ├── 00_fetch_ema_ma_dates.R       # download + clean EMA dates
│   ├── 01_data_cleaning.R            # data prep and merging
│   ├── 02_time_to_access.R           # survival analysis
│   └── 03_concordance.R              # concordance + mixed-effects regression
├── figures/
│   ├── km_by_country.png
│   ├── cox_forest_plot.png
│   ├── time_distribution_by_country.png
│   ├── approval_rate_by_technology.png
│   ├── concordance_summary.png
│   ├── concordance_heatmap.png
│   └── glmer_forest_plot.png
├── outputs/
│   ├── 02_results.xlsx               # survival analysis numerical results (7 sheets)
│   └── 03_results.xlsx               # concordance numerical results (5 sheets)
└── README.md
```

---

## Reproducing the analysis

```r
# Install required packages (once)
install.packages(c(
  "survival", "survminer", "dplyr", "tidyr", "readr", "ggplot2",
  "here", "broom", "broom.mixed", "scales", "forcats", "stringr",
  "lubridate", "readxl", "lme4", "openxlsx"
))

# Run in order
source("notebooks/00_fetch_ema_ma_dates.R")   # downloads ~2 MB from EMA website
source("notebooks/01_data_cleaning.R")
source("notebooks/02_time_to_access.R")
source("notebooks/03_concordance.R")
```

Each script saves its figures to `figures/` (overwriting previous versions) and numerical results to `outputs/`.

---

## Connection to thesis work

This project is exploratory groundwork for a master's thesis on Nordic orphan drug market access at Karolinska Institutet (HEPM programme). The main thesis-relevant findings are:

- **Country of evaluation is a stronger predictor of time-to-access than technology type** — Denmark approves drugs ~2.85× faster than Norway, after adjusting for what the drug is.
- **49% cross-country discordance** is high enough to raise serious questions about HTA harmonisation in the Nordics, relevant to the EU HTA Regulation (Regulation (EU) 2021/2282) which mandates joint clinical assessments for orphan drugs from 2025.
- The **real-world evidence paradox** (faster access despite weaker evidence) is consistent with published literature on value-based access for ultra-rare diseases, where RCTs are often impossible and ethical considerations override standard cost-effectiveness thresholds.

---

*Data manually curated by Margarida Bettencourt (2026). All source documents are publicly available from TLV, Medicinrådet, and Nye Metoder websites.*
