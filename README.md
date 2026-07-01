# Nordic Orphan Drug Market Access

Survival analysis and cross-country concordance of orphan drug HTA decisions across Sweden, Denmark, and Norway.

**Dataset:** 51 drugs × 3 countries = 174 decision records, manually curated from public TLV, Medicinrådet, and Nye Metoder assessments (2018–2026).

---

## Questions

1. **Time to access** — after EMA approval, how long does each country take to reach a positive decision? Does this vary by severity tier or technology type?
2. **Cross-country concordance** — when the same drug is evaluated in all three countries, do they agree? What predicts a split decision?

---

## Methods

- Kaplan-Meier survival curves + log-rank test (`survival`)
- Cox proportional hazards regression
- Mixed-effects logistic regression with drug-level random intercept (`lme4`)
- All analyses in R

---

## Structure

```
data/raw/          ← source Excel (not modified)
notebooks/
  01_data_cleaning.R    ← load, fix types, merge EMA dates, export CSV
  02_time_to_access.R   ← Kaplan-Meier + Cox regression
  03_concordance.R      ← cross-country agreement + mixed-effects logistic
outputs/figures/   ← generated plots (not tracked in git)
```

---

## Data notes

- Decision dates available for ~136/174 records; remainder have year only
- ICER data available for Sweden only (TLV)
- EMA approval dates merged from EMA COMP Excel to compute time-to-decision intervals
- All source data is public; PDFs archived separately

---

## Status

Work in progress.
