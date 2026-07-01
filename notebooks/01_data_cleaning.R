# =============================================================================
# 01_data_cleaning.R
# Nordic Orphan Drug Market Access — Data Cleaning & Preparation
# =============================================================================
# Goal: Load the master dataset, fix types, patch known missing values,
#       merge EMA approval dates, compute time-to-decision, save clean CSV.
#
# Missing value policy: NA stays NA. No imputation, no zero-filling.
# Records with missing values are dropped only within the specific analysis
# that requires that variable — not removed globally.
#
# Input:  data/raw/nordic_orphan_master.xlsx  (Master Data sheet)
#         data/raw/ema_orphan_designations_raw.xlsx  (EMA COMP export)
# Output: data/nordic_orphan_clean.csv
# =============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(here)

# -----------------------------------------------------------------------------
# 1. Load raw HTA data
# -----------------------------------------------------------------------------
raw <- read_excel(
  here("data", "raw", "nordic_orphan_master.xlsx"),
  sheet = "Master Data"
)

cat("Raw data dimensions:", nrow(raw), "rows x", ncol(raw), "columns\n")
cat("Countries:", paste(unique(raw$Country), collapse = ", "), "\n")
cat("Drugs:", n_distinct(raw$Drug), "\n\n")

# -----------------------------------------------------------------------------
# 2. Fix column types
# -----------------------------------------------------------------------------
df <- raw |>
  mutate(
    Is_Positive   = as.logical(Is_Positive),
    Is_Assessed   = as.logical(Is_Assessed),
    Severity_Tier = factor(Severity_Tier, levels = c(
      "Ultra-severe (immediately life-threatening)",
      "Severe (significantly life-limiting)",
      "Moderate (serious but not immediately life-threatening)"
    ), ordered = TRUE),
    Country             = factor(Country),
    Technology_Group    = factor(Technology_Group),
    Evidence_Type       = factor(Evidence_Type),
    Restriction_Level   = factor(Restriction_Level),
    Comparator_Category = factor(Comparator_Category),
    Decision_Year       = as.integer(Decision_Year)
  )

# -----------------------------------------------------------------------------
# 3. Patch known missing Decision_Year values
#    Source: verified from official HTA body websites (public records only)
# -----------------------------------------------------------------------------
df <- df |>
  mutate(Decision_Year = case_when(
    Drug == "Brineura" & Country == "Sweden" & is.na(Decision_Year) ~ 2018L,
    Drug == "Enspryng" & Country == "Norway" & is.na(Decision_Year) ~ 2021L,
    Drug == "Zynteglo" & Country == "Norway" & is.na(Decision_Year) ~ 2020L,
    Drug == "Lojuxta"  & Country == "Sweden" & is.na(Decision_Year) ~ 2025L,
    TRUE ~ Decision_Year
  ))

cat("Decision_Year patch applied. Missing after patch:",
    sum(is.na(df$Decision_Year)), "(was 38)\n\n")

# -----------------------------------------------------------------------------
# 4. Parse HTA decision date
#    Use Decision_Date_Raw where it looks like a real date.
#    For records with only Decision_Year, use July 1 as mid-year estimate.
#    This approximation is documented — do not over-interpret month-level
#    precision in the survival analysis.
# -----------------------------------------------------------------------------
df <- df |>
  mutate(
    # Try to parse the raw date string
    hta_date_parsed = suppressWarnings(as_date(Decision_Date_Raw)),

    # Where parsing failed but we have a year, use July 1 of that year
    hta_decision_date = case_when(
      !is.na(hta_date_parsed)              ~ hta_date_parsed,
      !is.na(Decision_Year)                ~ as_date(paste0(Decision_Year, "-07-01")),
      TRUE                                 ~ NA_Date_
    )
  ) |>
  select(-hta_date_parsed)

cat("HTA decision dates resolved:",
    sum(!is.na(df$hta_decision_date)), "/ 174\n\n")

# -----------------------------------------------------------------------------
# 5. Merge EMA first marketing authorisation dates
#    Source: ema_ma_dates.csv, built by 00_fetch_ema_ma_dates.R
#    This uses the "Marketing authorisation date" column from the EMA medicines
#    data table — the true first-authorisation date, not the COMP designation
#    date or a recent amendment date.
#    Join key: INN (cleaned to lowercase, trimmed)
# -----------------------------------------------------------------------------
ema_dates <- read.csv(
  here("data", "raw", "ema_ma_dates.csv"),
  stringsAsFactors = FALSE
) |>
  mutate(ema_ma_date = as_date(ema_ma_date)) |>
  select(inn_clean, ema_ma_date)

cat("EMA MA dates loaded:", nrow(ema_dates), "unique INNs\n")

# Join onto master data
df <- df |>
  mutate(inn_clean = str_to_lower(str_trim(INN))) |>
  left_join(ema_dates, by = "inn_clean") |>
  select(-inn_clean)

cat("EMA MA date matched for:",
    sum(!is.na(df$ema_ma_date)), "/ 174 records\n\n")

# -----------------------------------------------------------------------------
# 6. Compute time-to-decision (years from EMA first MA to HTA decision)
# -----------------------------------------------------------------------------
df <- df |>
  mutate(
    time_to_decision_years = as.numeric(
      interval(ema_ma_date, hta_decision_date) / dyears(1)
    ),
    # Flag records where time is implausible (data quality check)
    time_flag = case_when(
      is.na(time_to_decision_years)   ~ "missing",
      time_to_decision_years < 0      ~ "negative (HTA before MA — check dates)",
      time_to_decision_years > 10     ~ "long (>10 years)",
      TRUE                            ~ "ok"
    )
  )

cat("Time-to-decision summary:\n")
print(table(df$time_flag))
cat("\n")
cat("Median time to decision (years):",
    round(median(df$time_to_decision_years, na.rm = TRUE), 2), "\n\n")

# -----------------------------------------------------------------------------
# 7. Document remaining missing values
# -----------------------------------------------------------------------------
missing_summary <- df |>
  summarise(across(everything(), \(x) sum(is.na(x)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") |>
  filter(n_missing > 0) |>
  arrange(desc(n_missing))

cat("Remaining missing values by column:\n")
print(missing_summary)
cat("\n")

# -----------------------------------------------------------------------------
# 8. Save clean dataset
# -----------------------------------------------------------------------------
write.csv(df, here("data", "nordic_orphan_clean.csv"), row.names = FALSE, na = "")

cat("Clean data saved to data/nordic_orphan_clean.csv\n")
cat("Rows:", nrow(df), "| Columns:", ncol(df), "\n")
