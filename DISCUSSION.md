# Discussion: Assumptions, Findings, Limitations & Next Steps

This document explains what the analysis does, what it assumes, what it found, and what to do next. It is written for a general reader — no prior statistics knowledge required. Every assumption is flagged; every methodological decision is explained.

---

## 1. What We Set Out to Understand

This project investigates **orphan drug market access in Nordic countries**. The central questions are:

1. After a drug receives EMA marketing authorisation, how long does it take for Nordic HTA bodies (TLV, Medicinrådet, SBU in Sweden; TLV/NT-rådet in Denmark; Nye Metoder in Norway) to make a **positive** reimbursement decision?
2. Do the three countries differ significantly in how fast and how often they reimburse orphan drugs?
3. When the same drug is assessed in all three countries, do they tend to **agree**?
4. What drug or evidence characteristics predict a positive outcome?

This is important because: orphan drugs treat rare, severe diseases; patients in different Nordic countries can face very different delays or refusals for the same drug; and the decision-making factors are not always publicly explained.

---

## 2. Data Sources and Construction

| Source | What It Provides | Limitation |
|--------|-----------------|------------|
| Manual collection from TLV, NT-rådet, Nye Metoder public databases | Drug names, INN, decision outcomes, decision dates, technology type, evidence type, severity tier | Manual work = possible entry errors; some fields not systematically published |
| EMA medicines authorisation table (downloaded programmatically) | First EMA marketing authorisation date per drug | EMA file uses serial date numbers for some entries; solved with `origin = "1899-12-30"` |
| EMA COMP table (orphan designation) | NOT used for dates — only identifies which drugs had orphan designation | COMP file contains amendment dates, not first MA date — would produce incorrect time calculations |

The dataset covers **51 unique orphan drugs**, each assessed in one to three Nordic countries, giving **174 drug-country assessment records**.

---

## 3. Assumptions, Flagged

Each assumption is numbered and referenced in the R scripts.

### [A1] July 1 date approximation
When only a year is available for the HTA decision date, we use `July 1` as the mid-year estimate. This affects approximately 34 records.

- **Why**: 34 records had only a year in the source data, not an exact date.
- **Impact**: Adds up to ±6 months of noise to `time_to_decision_years`. With wide confidence intervals throughout (see Section 5), this is unlikely to change any conclusions.
- **Disclosure**: Approximately 20% of records are date-approximated. Any time-based finding should be interpreted with this caveat in mind.

### [A2] EMA first MA date as "time zero"
We measure time from EMA marketing authorisation to HTA decision, not from orphan designation or the date of submission to the HTA body.

- **Why**: EMA MA date is the earliest legally possible moment for a drug to be submitted for reimbursement. It is publicly available and reproducible. Submission dates to each HTA body are not systematically published.
- **Alternative considered**: Orphan designation date (from COMP file). Rejected because COMP date reflects when Eur Commission granted research incentives — which can precede the actual drug approval by years. Using it would grossly overstate the "waiting" time.
- **Implication**: Our time-to-access measures "how long after European approval did Nordic patients have to wait." It does not measure "how long after national submission."

### [A3] Negative decisions treated as censored, not as failures
Drugs that received a negative or restricted HTA decision are treated as **censored** in the survival analysis, not as permanent failures.

- **Why**: A censored observation means "we don't know the final outcome — the event might still happen." Drugs can be resubmitted after a negative decision, and sometimes get approved on the second attempt (Evrysdi in Norway is an example in this dataset). Treating negative decisions as permanent failures would ignore this possibility.
- **Consequence**: This slightly **underestimates** the true censoring time for drugs that were definitively rejected and never resubmitted. This is accepted as a conservative simplification.

### [A4] Concordance = current/final status
For cross-country concordance, we define each drug's status per country as the **best outcome ever achieved** (`max(Is_Positive)`). If a drug was rejected in 2021 and approved in 2023, it counts as positive.

- **Why**: From a patient access perspective, current reimbursement status is what matters. A drug that eventually got approved IS accessible.
- **Alternative**: Use the FIRST decision per drug-country pair (could show early rejection + later approval patterns). Not implemented here but relevant for studying re-submission dynamics.

### [A5] Clustered standard errors in Cox regression
The same drug appears in up to 3 countries → observations are correlated. We use `cluster = Drug` in `coxph()` to produce robust standard errors.

- **Why**: Standard Cox assumes independence between observations. Here, Spinraza in Norway and Spinraza in Sweden share the same clinical evidence package, manufacturer, and disease indication. Ignoring this would make SEs artificially small → overconfident p-values.
- **Technical note**: `tidy()` from broom does NOT propagate clustered CIs — it uses naive SEs. This is a known limitation. Point estimates (HRs) are correct; CIs are slightly too narrow.

### [A6] Exclusion of "Other biologic" and Evidence_Type from glmer models
Two predictors were excluded from the logistic regression models due to **quasi-complete separation**.

- **What is quasi-complete separation?** When a predictor perfectly (or near-perfectly) predicts the outcome, the maximum likelihood estimate is infinite — the model "learns" that this category is a guaranteed predictor. With small N, this happens by chance: 3 "Other biologic" drugs all happened to be approved, but that doesn't mean Other biologics always get approved — we just don't have negative examples.
- **Other biologic**: Only 2–3 records, all positive → excluded from regression. Reported in descriptive tables.
- **Evidence_Type "Real-world/registry"**: All drugs in this evidence category were approved in both Model A and Model B subsets → excluded from both models. Reported descriptively.
- **This is standard practice** in clinical and health economics statistics for small N. The key is documenting the exclusion (as here), not silently hiding it.

### [A7] Random intercept per drug in glmer
The `(1 | Drug)` term absorbs unmeasured drug-level heterogeneity.

- **What this means**: Everything specific to a drug that isn't captured by our covariates (exact efficacy magnitude, price negotiations, manufacturer's willingness to negotiate, clinical comparator strength, disease severity within tiers) is absorbed into the drug-specific random effect. This makes the fixed effects (Country, Technology_Group) cleaner estimates of those specific contributions.
- **Assumption**: Random effects are normally distributed with mean 0. This is untestable with 51 drugs but is the standard mixed-effects assumption.

---

## 4. Key Findings

### Time to access
- Median time from EMA authorisation to a positive HTA decision ranges from approximately **1.5 to 3.0 years** across the three countries. (Exact values in `outputs/02_results.xlsx`, sheet "Country Summary".)
- The **log-rank test** comparing the three countries' KM curves provides a p-value for whether the time-to-access distributions are statistically different. Given N~129, the test is likely underpowered to detect moderate differences — interpret directionally, not definitively.

### Approval rates
- Approval rates differ across countries and across technology types. **Biologics** and **cell therapies** tend to have higher approval rates than gene therapies, consistent with established clinical evidence bases and pricing frameworks for the former.
- Country differences in approval rate reflect not just clinical judgement but also **health economic frameworks**: Sweden's TLV applies a cost-per-QALY threshold; Denmark's Medicinrådet uses managed entry; Norway's Nye Metoder has explicit budget impact considerations.

### Concordance
- A meaningful proportion of drugs evaluated in all three countries receive **discordant decisions** — one country approves while another rejects.
- Discordance tends to concentrate in drugs with **higher uncertainty** in clinical evidence (small trials, surrogate endpoints) and/or **high cost** per patient.
- The heatmap (`figures/concordance_heatmap.png`) shows which specific drugs are discordant — useful for qualitative case-study follow-up.

### Regression results
- Country effects (Norway vs Denmark vs Sweden) and Technology_Group effects are estimable, but **CIs are wide**. No individual predictor is clearly dominant in either the Cox or glmer models.
- This is expected: with 51 drugs and several correlated predictors, we are underpowered for formal statistical discrimination. **The models are exploratory, not confirmatory.**
- The random effect variance in glmer (`Variance` under `Random effects` in the model summary) indicates how much between-drug variation exists AFTER accounting for country and technology. A large variance = drugs differ a lot in their base approval probability, net of covariates.

---

## 5. Limitations

**Sample size**: 51 drugs is small for multivariate regression. This is a consequence of the orphan drug space — by definition, few drugs. All statistical results should be treated as hypothesis-generating, not as definitive evidence.

**Selection bias**: We only included drugs with EMA orphan designation that WENT THROUGH HTA in at least one Nordic country. Drugs that never reached HTA (e.g. manufacturers never submitted due to anticipated rejection) are excluded. This understates the access problem.

**Missing Severity_Tier**: ~67% of records are missing this variable. It was coded manually from HTA body decision documents, which do not use a standardised severity classification. Model B (which includes Severity_Tier) covers only ~33% of the dataset.

**Time period**: The dataset covers a specific historical window. Regulatory frameworks, HTA methodologies, and pricing policies have changed over this period. Trends over time are not explicitly modelled.

**Censoring mechanism**: We assume censoring is non-informative — that whether a drug is censored is unrelated to when it would have received a positive decision. This may not hold if, for example, drugs with no chance of approval are more likely to be re-submitted late.

**Manual data entry**: The master dataset was compiled manually from public HTA databases. Despite validation checks in `01_data_cleaning.R`, some data entry errors may remain.

---

## 6. Lessons Learned (Technical)

**EMA date sourcing**: The COMP file (orphan designations) stores AMENDMENT dates, not first marketing authorisation dates. This took debugging to discover and means you CANNOT use COMP for time-zero calculations. Use the EMA medicines authorisation table instead (`Marketing authorisation date` column, header at row 9 → `skip = 8` in readxl).

**pivot_wider list-column error**: If a drug has multiple HTA records per country (rejection → re-submission → approval), `pivot_wider()` fails with "values are not uniquely identified." Always deduplicate to one row per Drug × Country before pivoting.

**droplevels() after filter**: Subsetting a data frame does NOT remove empty factor levels. If you filter to a subset and then fit a model, R will try to estimate coefficients for factor levels that no longer exist in the subset → NA or errors. Always call `droplevels()` after filtering.

**Quasi-complete separation detection**: If a glmer or glm model produces `estimate: 12.3, std.error: 143.6`, that is separation. Don't report those coefficients — exclude the offending predictor or merge its levels and document the decision.

**ggforest() incompatibility with cluster =**: The survminer `ggforest()` function does not support Cox models fitted with `cluster =`. Use `broom::tidy()` + manual ggplot2 instead (see `02_time_to_access.R`, Section 5a).

**ggsurvplot() requires png()/dev.off()**: `ggsurvplot()` objects contain two ggplot panels combined with `grid.arrange()` internally. `ggsave()` only handles single ggplot objects. To save a survminer plot to file, wrap it in `png()` / `print()` / `dev.off()`.

---

## 7. Next Steps

### Immediate (complete before thesis submission)
1. **Validate the master dataset** against a second reviewer or source document spot-check. Priority: the 5 records with negative `time_to_decision_years` (HTA before EMA MA) and the 4 manually patched `Decision_Year` values.
2. **Recover Severity_Tier for more records**: the current 33% coverage limits the Severity model substantially. Decision documents often mention severity indirectly even when not using standardised language.

### Extensions for thesis
3. **Time trend analysis**: Are orphan drugs being reimbursed faster or slower over time? Are earlier-approved drugs (pre-2015) systematically different from recent ones? Add a `cohort` variable splitting by EMA approval year and test for interaction.
4. **Restriction level analysis**: The dataset includes a `Restriction_Level` variable (unrestricted, hospital only, specific subpopulation). Analysing which drugs get conditions/restrictions (vs clean approvals) is a separate and policy-relevant question.
5. **Joint Nordic HTA**: Nine drugs have been evaluated through the FINOSE joint Nordic HTA process. Add a flag for joint vs national-only assessment and test whether joint HTA is associated with better concordance or faster access.

### Future projects (portfolio)
6. **Automated HTA scraper** (Python): Build a scraper for TLV, Nye Metoder, and Medicinrådet decision PDFs to automate future data collection. NLP-based classification of severity and restriction level from free text.
7. **Full EMA orphan drug screening** (Python): Cross-reference all EMA-designated orphan drugs against TLV/NT decisions to quantify what proportion of EMA-authorised orphan drugs have NEVER reached Nordic HTA — an uncharted shadow of the access problem.
8. **Swedish MEA descriptive analysis** (R): TLV publishes managed-entry agreement (MEA/risk-sharing) data. A "then vs now" analysis comparing MEA prevalence and structure across time periods would be a natural third portfolio project.

---

## 8. Reproducibility

To reproduce the full analysis from scratch:

```r
# In R, with working directory set to the repo root:
source("notebooks/00_fetch_ema_ma_dates.R")   # downloads EMA data, creates data/raw/ema_ma_dates.csv
source("notebooks/01_data_cleaning.R")          # cleans master dataset, creates data/nordic_orphan_clean.csv
source("notebooks/02_time_to_access.R")         # survival analysis, outputs to figures/ and outputs/
source("notebooks/03_concordance.R")             # concordance + regression, outputs to figures/ and outputs/
```

Required R packages: `readxl`, `dplyr`, `lubridate`, `stringr`, `tidyr`, `here`, `survival`, `survminer`, `ggplot2`, `broom`, `scales`, `openxlsx`, `lme4`, `broom.mixed`, `forcats`, `readr`, `fuzzyjoin`, `curl`.

Install all with:
```r
install.packages(c("readxl", "dplyr", "lubridate", "stringr", "tidyr", "here",
                   "survival", "survminer", "ggplot2", "broom", "scales",
                   "openxlsx", "lme4", "broom.mixed", "forcats", "readr",
                   "fuzzyjoin", "curl"))
```
