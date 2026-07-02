# =============================================================================
# 01_data_cleaning.R
# Load, clean, and prepare the Nordic orphan drug HTA dataset
# =============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Loads the raw master dataset from Excel, fixes column types, fills a few
# known missing values from public sources, merges EMA first-authorisation dates,
# computes time-to-HTA-decision in years, and saves a clean CSV for analysis.
#
# MISSING VALUE POLICY
# --------------------
# NA stays NA. We do not impute, zero-fill, or average-fill missing values.
# Records with missing values are dropped ONLY within the specific analysis
# that requires that variable — not removed from the dataset globally.
# Rationale: imputation would introduce bias in a small dataset where ~30% of
# records are missing Severity_Tier and ~22% are missing exact decision dates.
# Better to be honest about what we don't know.
#
# ASSUMPTIONS
# -----------
# [A1] Decision_Date_Raw is treated as a calendar date when it looks like one.
#      When only Decision_Year is available, we use July 1 as a mid-year estimate.
#      This adds ~±6 months of noise to the survival analysis — acceptable given
#      the wide confidence intervals that already exist at N=51 drugs.
# [A2] EMA marketing authorisation date is the correct "time 0" for time-to-access.
#      We use the first EMA MA date (from script 00), not the orphan designation
#      date, not the submission date, and not any individual country's approval date.
# [A3] time_to_decision_years < 0 means the HTA decision date precedes the EMA
#      MA date — impossible in theory. These 5 records likely reflect data entry
#      issues (wrong year in Decision_Year) or drugs that were assessed before the
#      final EMA MA was formalised. They are flagged "negative" and excluded from
#      survival analysis. We do NOT delete them — they stay in the clean dataset
#      for transparency.
# [A4] Records with time_to_decision_years > 10 are flagged "long" but kept.
#      Some older drugs (e.g. Fabrazyme, approved by EMA in 2001) genuinely did
#      have late Nordic HTA assessments. This is a real pattern, not an error.
#
# Input:  data/raw/nordic_orphan_master.xlsx (sheet: "Master Data")
#         data/raw/ema_ma_dates.csv (from 00_fetch_ema_ma_dates.R)
# Output: data/nordic_orphan_clean.csv
# =============================================================================

library(readxl)    # read Excel files
library(dplyr)     # data manipulation
library(lubridate) # date arithmetic (interval, dyears, as_date)
library(stringr)   # string cleaning (str_to_lower, str_trim)
library(tidyr)     # pivot_longer (for missing value summary)
library(here)      # relative paths anchored at project root

# -----------------------------------------------------------------------------
# 1. Load raw HTA data
# -----------------------------------------------------------------------------
raw <- read_excel(
  here("data", "raw", "nordic_orphan_master.xlsx"),
  sheet = "Master Data"
)

cat("Raw data loaded:", nrow(raw), "rows x", ncol(raw), "columns\n")
cat("Countries:", paste(unique(raw$Country), collapse = ", "), "\n")
cat("Unique drugs:", n_distinct(raw$Drug), "\n\n")

# -----------------------------------------------------------------------------
# 2. Fix column types
# -----------------------------------------------------------------------------
# Excel doesn't distinguish TRUE/FALSE from 1/0 or "Yes"/"No".
# R reads everything as character or numeric. We convert explicitly.
#
# Severity_Tier is an ORDERED factor: Ultra-severe > Severe > Moderate.
# This ordering matters for any analysis that treats severity as numeric,
# e.g. ordinal regression or rank-based tests.

df <- raw |>
  mutate(
    # Binary flags
    Is_Positive = as.logical(Is_Positive),
    Is_Assessed = as.logical(Is_Assessed),

    # Ordered factor: severity from highest to lowest need
    # (Ultra-severe = immediately life-threatening, no treatment alternatives)
    Severity_Tier = factor(Severity_Tier, levels = c(
      "Ultra-severe (immediately life-threatening)",
      "Severe (significantly life-limiting)",
      "Moderate (serious but not immediately life-threatening)"
    ), ordered = TRUE),

    # Unordered factors — no inherent ordering assumed
    Country             = factor(Country),
    Technology_Group    = factor(Technology_Group),
    Evidence_Type       = factor(Evidence_Type),
    Restriction_Level   = factor(Restriction_Level),
    Comparator_Category = factor(Comparator_Category),

    Decision_Year = as.integer(Decision_Year)
  )

# -----------------------------------------------------------------------------
# 3. Patch known missing Decision_Year values
# -----------------------------------------------------------------------------
# A handful of records have no decision year in the master data.
# We fill these from the original HTA body decision documents (public sources).
# Each patch is documented below with the drug name and country.
# NOTE: these are hard-coded fixes — if the master data is updated, check these.

df <- df |>
  mutate(Decision_Year = case_when(
    Drug == "Brineura" & Country == "Sweden"  & is.na(Decision_Year) ~ 2018L,  # TLV decision
    Drug == "Enspryng" & Country == "Norway"  & is.na(Decision_Year) ~ 2021L,  # Nye Metoder
    Drug == "Zynteglo" & Country == "Norway"  & is.na(Decision_Year) ~ 2020L,  # Nye Metoder
    Drug == "Lojuxta"  & Country == "Sweden"  & is.na(Decision_Year) ~ 2025L,  # TLV subvention
    TRUE ~ Decision_Year
  ))

cat("Decision_Year patches applied. Missing Decision_Year after patch:",
    sum(is.na(df$Decision_Year)), "\n\n")

# -----------------------------------------------------------------------------
# 4. Parse HTA decision date
# -----------------------------------------------------------------------------
# ASSUMPTION [A1]: July 1 approximation for records where we only have the year.
#
# The raw date column (Decision_Date_Raw) contains a mix of:
#   - Proper ISO dates: "2021-02-24"
#   - Descriptive strings: "Original recommendation: 26 May 2021; follow-up..."
#   - Just a year: "2021"
#
# as_date() from lubridate parses well-formatted dates and returns NA for
# anything it can't parse. suppressWarnings() silences the "failed to parse"
# messages for the non-date strings — we handle those via Decision_Year.
#
# The July 1 approximation introduces at most ±6 months of error on the
# time_to_decision_years variable. Given the wide confidence intervals in all
# our models (N is small), this is unlikely to change any conclusions.

df <- df |>
  mutate(
    # Try to parse the raw date string — will be NA for descriptive text
    hta_date_parsed = suppressWarnings(as_date(Decision_Date_Raw)),

    # Use parsed date if available; otherwise July 1 of the year
    hta_decision_date = case_when(
      !is.na(hta_date_parsed) ~ hta_date_parsed,
      !is.na(Decision_Year)   ~ as_date(paste0(Decision_Year, "-07-01")),
      TRUE                    ~ NA_Date_
    )
  ) |>
  select(-hta_date_parsed)  # remove intermediate column

cat("HTA decision dates resolved:", sum(!is.na(df$hta_decision_date)), "/ 174\n")
cat("  (dates approximated as July 1 for", sum(!is.na(df$Decision_Year) &
    is.na(suppressWarnings(as_date(df$Decision_Date_Raw)))), "records)\n\n")

# -----------------------------------------------------------------------------
# 5. Merge EMA first marketing authorisation dates
# -----------------------------------------------------------------------------
# ASSUMPTION [A2]: EMA MA date = time 0 for our time-to-decision calculation.
#
# WHY NOT use the EMA COMP orphan designation date?
# Because COMP dates track when a drug received ORPHAN DESIGNATION (a research
# incentive, not marketing approval). Drugs can have orphan designation for years
# before they're authorised for sale. Using COMP dates would massively overstate
# the time "waiting" for HTA decisions.
#
# WHY NOT use each country's submission date?
# Because submission dates are not systematically publicly available for all
# drugs and countries. Using EMA MA date is a principled, reproducible choice —
# it represents when the product was available to be submitted for reimbursement.
#
# The join is on inn_clean (INN standardised to lowercase, trimmed).

ema_dates <- read.csv(
  here("data", "raw", "ema_ma_dates.csv"),
  stringsAsFactors = FALSE
) |>
  mutate(ema_ma_date = as_date(ema_ma_date)) |>
  select(inn_clean, ema_ma_date)

cat("EMA MA dates loaded:", nrow(ema_dates), "unique INNs\n")

df <- df |>
  mutate(inn_clean = str_to_lower(str_trim(INN))) |>
  left_join(ema_dates, by = "inn_clean") |>
  select(-inn_clean)  # remove the temporary join key

cat("EMA MA date matched for:", sum(!is.na(df$ema_ma_date)), "/ 174 records\n")
cat("Unmatched drugs (EMA date will be NA):\n")
df |> filter(is.na(ema_ma_date)) |> distinct(Drug) |> pull(Drug) |> print()
cat("\n")

# -----------------------------------------------------------------------------
# 6. Compute time-to-decision (years from EMA first MA to HTA decision)
# -----------------------------------------------------------------------------
# ASSUMPTION [A3 & A4]: See script header.
#
# We use lubridate::interval() / dyears() instead of simple subtraction because:
# - dyears(1) = 365.25 days (accounts for leap years)
# - interval() correctly handles month-length differences
# The result is a decimal — e.g. 2.03 years = just over 2 years.
#
# The time_flag column is a data-quality audit trail:
# "ok"       = plausible (0–10 years)
# "missing"  = either EMA date or HTA date is NA
# "negative" = HTA decision appears to predate EMA authorisation — data issue
# "long"     = >10 years — unusual but possible for older drugs

df <- df |>
  mutate(
    time_to_decision_years = as.numeric(
      interval(ema_ma_date, hta_decision_date) / dyears(1)
    ),
    time_flag = case_when(
      is.na(time_to_decision_years)  ~ "missing",
      time_to_decision_years < 0     ~ "negative (HTA before MA — check dates)",
      time_to_decision_years > 10    ~ "long (>10 years)",
      TRUE                           ~ "ok"
    )
  )

cat("Time-to-decision flag summary:\n")
print(table(df$time_flag))
cat("\nMedian time (records flagged 'ok' only):",
    round(median(df$time_to_decision_years[df$time_flag == "ok"], na.rm = TRUE), 2),
    "years\n\n")

# -----------------------------------------------------------------------------
# 7. Missing value audit
# -----------------------------------------------------------------------------
# Print a full accounting of what's missing and why.
# This is the foundation for the analysis-level filter decisions:
# - Survival analysis (script 02) uses records with time_flag == "ok"
# - Concordance analysis (script 03) uses all records with Is_Positive not NA
# - Severity_Tier is missing for ~70% — only included where available

missing_summary <- df |>
  summarise(across(everything(), \(x) sum(is.na(x)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") |>
  filter(n_missing > 0) |>
  arrange(desc(n_missing))

cat("Remaining missing values by column:\n")
print(missing_summary)
cat("\nNote: ICER columns are Sweden-only by design. Severity_Tier coded from\n")
cat("decision text — not systematically reported by all HTA bodies.\n\n")

# -----------------------------------------------------------------------------
# 8. Save clean dataset
# -----------------------------------------------------------------------------
write.csv(df, here("data", "nordic_orphan_clean.csv"), row.names = FALSE, na = "")

cat("Saved: data/nordic_orphan_clean.csv\n")
cat("Rows:", nrow(df), "| Columns:", ncol(df), "\n")
cat("Ready for analysis in scripts 02 and 03.\n")
