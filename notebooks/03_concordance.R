# =============================================================================
# 03_concordance.R
# Nordic Orphan Drug Market Access — Cross-Country Concordance Analysis
# =============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Asks: when the same orphan drug is evaluated in TWO OR THREE Nordic countries,
# do they reach the same reimbursement decision? And what predicts a positive
# decision in the first place?
#
# TWO COMPLEMENTARY ANALYSES
# --------------------------
# (A) DESCRIPTIVE CONCORDANCE
#     Reshapes the data from "one row per drug-country assessment" to
#     "one row per drug", then classifies each drug as:
#       - Concordant: all countries that assessed it agreed
#       - Discordant: at least one country said yes and another said no
#       - Single country only: no comparison possible
#     Also visualises which drugs are discordant in a heatmap.
#
# (B) MIXED-EFFECTS LOGISTIC REGRESSION (glmer)
#     Models what DRIVES positive decisions across all records (back to the
#     long format, one row per drug-country). Logistic regression estimates
#     how Country, Technology_Group, etc. affect the odds of approval.
#     Mixed-effects accounts for the fact that the same drug appears in multiple
#     countries (within-drug correlation).
#
# WHY MIXED-EFFECTS, NOT PLAIN LOGISTIC REGRESSION?
# --------------------------------------------------
# In plain logistic regression, observations are assumed independent.
# But Spinraza in Norway and Spinraza in Sweden are NOT independent — they
# share the same underlying drug characteristics (clinical evidence, mechanism,
# disease indication, manufacturer). If we ignore this, we:
# (a) double/triple-count the drug's contribution to the parameter estimates
# (b) get SEs that are too small → false precision
#
# Mixed-effects logistic regression adds a random intercept per drug:
#   Is_Positive ~ Country + Technology_Group + ... + (1 | Drug)
#
# The (1 | Drug) term absorbs all unmeasured drug-level factors — clinical benefit,
# price, how robust the evidence package is, disease severity, manufacturer
# strategy — and allows the FIXED effects (Country, Technology) to represent
# country-level and technology-level contributions cleanly, net of those drug
# factors.
#
# ASSUMPTIONS
# -----------
# [A1] For concordance, we take the CURRENT status of each drug per country.
#      If a drug was rejected in 2021 and then approved in 2024, we count it
#      as POSITIVE. This uses max(Is_Positive) — see Section 2 below.
# [A2] The mixed-effects model assumes the drug random effects are normally
#      distributed with mean 0. This is a standard assumption in lme4.
# [A3] Wald confidence intervals for glmer odds ratios (from broom.mixed)
#      are approximate. Exact CIs from likelihood profiling exist but are
#      computationally intensive and overkill for N~170 exploratory analysis.
# [A4] Quasi-complete separation: with small N, some predictor categories
#      perfectly predict the outcome → infinite coefficient. We exclude
#      "Other biologic" and "Evidence_Type" from the regression. See Section 5.
#
# Input:  data/nordic_orphan_clean.csv (from 01_data_cleaning.R)
# Output: figures/concordance_summary.png  — bar chart: concordance overview
#         figures/concordance_heatmap.png  — drug × country outcome grid
#         figures/glmer_forest_plot.png    — odds ratios from mixed-effects model
#         outputs/03_results.xlsx          — all numerical results
# =============================================================================

# Load packages ---------------------------------------------------------------
library(dplyr)       # data manipulation
library(tidyr)       # pivot_wider (reshape long → wide per drug)
library(readr)       # read_csv
library(ggplot2)     # visualisation
library(here)        # project-relative paths
library(lme4)        # glmer() — mixed-effects logistic regression
library(broom.mixed) # tidy() for glmer objects (base broom doesn't handle mixed models)
library(scales)      # percent_format() on plot axes
library(forcats)     # fct_rev() to reverse factor ordering in plots
library(stringr)     # str_remove() for cleaning label names
library(openxlsx)    # Excel export

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)

# =============================================================================
# SECTION 1: Load data
# =============================================================================
df <- read_csv(here("data", "nordic_orphan_clean.csv"), show_col_types = FALSE) |>
  mutate(
    # Encode variables as factors.
    # Norway = reference for Country (most data, most commonly benchmarked)
    Country          = factor(Country, levels = c("Norway", "Denmark", "Sweden")),
    Technology_Group = factor(Technology_Group),
    Evidence_Type    = factor(Evidence_Type),

    # Convert logical Is_Positive to integer (1 = positive, 0 = negative/NA)
    # Required for glmer() which expects a numeric response in binomial models
    Is_Positive      = as.integer(Is_Positive)
  )

cat("Total records loaded:", nrow(df), "\n")
cat("Drugs:", n_distinct(df$Drug), "\n\n")

# =============================================================================
# SECTION 2: Reshape to wide format — one row per drug
# =============================================================================
#
# WHY DEDUPLICATION FIRST?
# Some drugs have multiple HTA records per country:
# e.g. Drug X rejected in Norway (2020), then resubmitted and approved (2023).
# In the raw data, this is two rows for Norway. If we pivot to wide without
# deduplicating, we'd get a list-column (two values for the same country column)
# which crashes pivot_wider with "values are not uniquely identified".
#
# FIX: for concordance purposes, we want the CURRENT status of each drug in
# each country. We define this as: if the drug was EVER positively reimbursed,
# current status = 1. This uses max(Is_Positive) per Drug × Country pair.
# max() treats 1 > 0, so if any record is positive, max = 1.
#
# ASSUMPTION [A1]: "Final" status = best outcome ever achieved, not latest.
# Rationale: from a patient access perspective, a drug that was rejected in
# 2020 then approved in 2023 IS accessible today. We want to capture that.

# Step 1: extract drug-level features (same for all countries)
# slice(1) takes the first row per drug — since these columns are the same
# across all country assessments for a given drug, the choice of row doesn't matter
drug_features <- df |>
  group_by(Drug) |>
  slice(1) |>
  ungroup() |>
  select(Drug, Technology_Group, Severity_Tier, Therapeutic_Area)

# Step 2: deduplicate — one row per Drug × Country, taking max Is_Positive
df_dedup <- df |>
  filter(!is.na(Is_Positive)) |>             # only records with a known outcome
  group_by(Drug, Country) |>
  summarise(Is_Positive = max(Is_Positive), .groups = "drop")  # current status

cat("Records after deduplication (one per Drug × Country):", nrow(df_dedup), "\n\n")

# Step 3: pivot from long (rows = drug-country) to wide (rows = drug, cols = countries)
# names_from = Country → creates column for each country
# values_from = Is_Positive → fills each country column with the outcome
# names_prefix = "pos_" → avoids column names starting with "Norway" (bad R syntax)
df_wide <- df_dedup |>
  pivot_wider(
    names_from   = Country,
    values_from  = Is_Positive,
    names_prefix = "pos_"
  ) |>
  left_join(drug_features, by = "Drug") |>    # add drug-level features back in
  mutate(
    # Count how many countries actually assessed this drug
    # across(starts_with("pos_")) selects all 3 country columns at once
    n_countries  = rowSums(!is.na(across(starts_with("pos_")))),

    # Count how many of those assessments were positive
    # na.rm = TRUE: countries where drug wasn't assessed (NA) don't count as negative
    n_positive   = rowSums(across(starts_with("pos_")), na.rm = TRUE),

    # Classify concordance:
    # "Concordant" = all countries that assessed it AGREED (all yes or all no)
    # "Discordant" = at least one yes AND at least one no
    # "Single country only" = can't assess concordance with only one observation
    concordance  = case_when(
      n_countries == 1 ~ "Single country only",
      n_countries == 2 & n_positive == 2 ~ "Concordant — both positive",
      n_countries == 2 & n_positive == 0 ~ "Concordant — both negative",
      n_countries == 2 & n_positive == 1 ~ "Discordant",
      n_countries == 3 & n_positive == 3 ~ "Concordant — all positive",
      n_countries == 3 & n_positive == 0 ~ "Concordant — all negative",
      n_countries == 3 & n_positive %in% c(1, 2) ~ "Discordant",
      TRUE ~ NA_character_
    ),
    # Boolean flag: TRUE if concordant (any type), FALSE if discordant
    # Used for quick counts below
    is_concordant = concordance != "Discordant" & concordance != "Single country only"
  )

cat("Drugs by number of countries evaluated:\n")
print(table(df_wide$n_countries))
cat("\nConcordance overview:\n")
print(table(df_wide$concordance, useNA = "ifany"))
cat("\n")

# Report on the all-3-country subset separately — this is the most informative
# subset because concordance is most meaningful when all three countries had a
# chance to weigh in
df_all3 <- df_wide |> filter(n_countries == 3)
cat("Drugs evaluated in all 3 countries:", nrow(df_all3), "\n")
cat("  Fully concordant:", sum(df_all3$is_concordant), "\n")
cat("  Discordant:      ", sum(!df_all3$is_concordant), "\n\n")

# Print the list of discordant drugs — useful for qualitative follow-up
# (e.g. which drug was approved in 2 countries but rejected in 1? Why?)
discordant_drugs <- df_wide |>
  filter(concordance == "Discordant") |>
  arrange(n_countries, Drug) |>
  select(Drug, Technology_Group, Severity_Tier,
         pos_Norway, pos_Denmark, pos_Sweden,
         n_countries, n_positive, concordance)

cat("Discordant drugs:\n")
print(discordant_drugs)
cat("\n")

# =============================================================================
# SECTION 3: Concordance summary bar chart
# =============================================================================
#
# Shows the breakdown of drugs into concordance categories, grouped by how
# many countries assessed the drug. Stacked horizontal bar = easy to compare
# counts. White numbers inside bars show the count for each category.
#
# COLOUR LOGIC:
# Green = positive concordance (all approved)
# Red = negative concordance (all rejected)
# Orange = discordant
# Grey = single country only (no comparison possible)
# Lighter shades = 2-country subgroups

concordance_levels <- c(
  "Concordant — all positive",
  "Concordant — both positive",
  "Concordant — all negative",
  "Concordant — both negative",
  "Discordant",
  "Single country only"
)
concordance_colours <- c(
  "Concordant — all positive"  = "#3DAA6D",
  "Concordant — both positive" = "#82C9A0",
  "Concordant — all negative"  = "#cc4444",
  "Concordant — both negative" = "#e89090",
  "Discordant"                 = "#E8A838",
  "Single country only"        = "#cccccc"
)

p_concordance <- df_wide |>
  mutate(
    concordance   = factor(concordance, levels = concordance_levels),
    # Human-readable label for the facet / y axis
    n_label       = paste0(n_countries, " countr", ifelse(n_countries == 1, "y", "ies"))
  ) |>
  count(n_label, concordance) |>
  ggplot(aes(x = fct_rev(n_label), y = n, fill = concordance)) +
  geom_col(width = 0.6) +
  # White bold labels inside each segment — use vjust = 0.5 for vertical centre
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = concordance_colours, drop = FALSE) +
  coord_flip() +    # horizontal bars are easier to read with long category names
  labs(
    title    = "Cross-country concordance of HTA decisions",
    subtitle = paste0(
      "N=", nrow(df_wide), " drugs; ",
      nrow(df_all3), " evaluated in all 3 countries, ",
      sum(df_all3$is_concordant), " fully concordant"
    ),
    x    = NULL,
    y    = "Number of drugs",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9)) +
  guides(fill = guide_legend(nrow = 2))   # legend in 2 rows to avoid overflow

print(p_concordance)
ggsave(here("figures", "concordance_summary.png"),
       plot = p_concordance, width = 9, height = 5, dpi = 150)

# =============================================================================
# SECTION 4: Outcome heatmap — all drugs × all countries
# =============================================================================
#
# A grid where rows = drugs, columns = countries, cells = outcome (green/red/grey).
# Sorted so fully concordant positives appear at the top, moving down through
# discordant and then concordant negatives.
#
# HOW TO READ THIS:
# - All three cells in a row the same colour = concordant drug
# - Mixed colours (green + red) = discordant drug
# - Grey = not assessed in that country (drug only went to 1 or 2 countries)
#
# WHY SORT BY CONCORDANCE, NOT ALPHABETICALLY?
# Alphabetical order is arbitrary — sorting by concordance category reveals the
# pattern structure: how many drugs are universally approved? How many are
# universally rejected? Where does disagreement cluster?

drug_order <- df_wide |>
  mutate(sort_key = case_when(
    concordance == "Concordant — all positive"  ~ 1,
    concordance == "Concordant — both positive" ~ 2,
    concordance == "Discordant"                 ~ 3,
    concordance == "Concordant — both negative" ~ 4,
    concordance == "Concordant — all negative"  ~ 5,
    TRUE                                        ~ 6   # single country only
  )) |>
  arrange(sort_key, Drug) |>   # within each concordance group, alphabetical
  pull(Drug)

# Convert Is_Positive to a labelled character for the fill scale
# (factor with explicit NA handling gives us control over the "not assessed" colour)
df_heatmap <- df |>
  select(Drug, Country, Is_Positive) |>
  mutate(
    outcome = case_when(
      Is_Positive == 1 ~ "Positive",
      Is_Positive == 0 ~ "Negative",
      TRUE             ~ NA_character_    # not assessed → grey
    ),
    Drug    = factor(Drug, levels = rev(drug_order)),   # rev() so top of plot = first drug
    Country = factor(Country, levels = c("Norway", "Denmark", "Sweden"))
  )

p_heatmap <- ggplot(df_heatmap, aes(x = Country, y = Drug, fill = outcome)) +
  geom_tile(colour = "white", linewidth = 0.5) +   # white lines between cells = grid
  scale_fill_manual(
    values   = c("Positive" = "#3DAA6D", "Negative" = "#cc4444"),
    na.value = "#e8e8e8",   # grey for not assessed
    name     = "HTA decision"
  ) +
  labs(
    title    = "HTA decision outcomes by drug and country",
    subtitle = "Sorted by concordance: all-positive → discordant → all-negative → not evaluated (grey)",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y      = element_text(size = 7),   # small text to fit 51 drugs
    legend.position  = "bottom",
    panel.grid       = element_blank()            # no grid lines — tiles ARE the grid
  )

print(p_heatmap)
# Tall figure: 51 drugs × enough height per row = 14 inches at 150dpi
ggsave(here("figures", "concordance_heatmap.png"),
       plot = p_heatmap, width = 7, height = 14, dpi = 150)

# =============================================================================
# SECTION 5: Mixed-effects logistic regression
# =============================================================================
#
# SETUP: LONG FORMAT (one row per drug-country assessment)
# The glmer model operates on the long-format data, NOT the wide format.
# Wide format was only needed for the concordance classification.

cat("=== Mixed-effects logistic regression ===\n\n")

df_glmer <- df |>
  filter(!is.na(Is_Positive)) |>
  mutate(
    Country          = factor(Country, levels = c("Norway", "Denmark", "Sweden")),
    Technology_Group = factor(Technology_Group),
    Evidence_Type    = factor(Evidence_Type),
    Drug             = factor(Drug)   # glmer needs Drug as factor for random effects
  )

# ---- WHY WE EXCLUDE "Other biologic" AND Evidence_Type ----------------------
#
# QUASI-COMPLETE SEPARATION: in logistic/mixed-effects regression, separation
# occurs when a predictor PERFECTLY predicts the outcome → the maximum likelihood
# estimate for that coefficient is +∞ or -∞, the SE is astronomically large,
# and the confidence intervals are meaningless (e.g. "OR = 1e9 [CI: 0, 1e15]").
#
# Two cases in our data:
#
# Case 1: "Other biologic" (Technology_Group)
#   This category contains only 2–3 drugs, and they all have the SAME outcome.
#   If all are positive: the model would estimate OR → +∞ ("being Other biologic
#   perfectly predicts approval"). But that's not a real finding — it's just that
#   our sample happened to have no Other biologic failures. With 2 drugs,
#   this is pure chance.
#   FIX: exclude "Other biologic" records. Report the raw n and outcomes in
#   descriptive tables so we're transparent about what's excluded.
#
# Case 2: Evidence_Type "Real-world / registry"
#   All drugs in this evidence category were approved → same problem as above.
#   In the Model A subset (N~171), all real-world evidence records are positive.
#   FIX: drop Evidence_Type from Model A entirely. The individual drugs and their
#   outcomes are still described in the heatmap and concordance tables.
#
# This kind of exclusion is STANDARD PRACTICE for small-N logistic models.
# The key is to document it (as we do here) rather than silently exclude or
# get nonsensical estimates.

other_bio_n <- sum(df_glmer$Technology_Group == "Other biologic", na.rm = TRUE)
cat("Note: 'Other biologic' has", other_bio_n,
    "records — excluded from regression (quasi-complete separation).\n")
cat("Note: Evidence_Type 'Real-world/registry' also causes separation — dropped from Model A.\n")
cat("      Both are reported in descriptive tables instead.\n\n")

# --- Model A: Country + Technology_Group + (1 | Drug), maximum N -------------
#
# After removing "Other biologic" records, we call droplevels() to tell R that
# the "Other biologic" factor level no longer exists in this subset.
# Without droplevels(), R would still try to estimate a coefficient for it
# (with 0 observations) → NA in the model output.

df_glmer_a <- df_glmer |>
  filter(Technology_Group != "Other biologic") |>
  mutate(Technology_Group = droplevels(Technology_Group))

cat("Model A — Country + Technology_Group + (1|Drug)\n")
cat("N =", nrow(df_glmer_a), "records,", n_distinct(df_glmer_a$Drug), "drugs\n\n")

# optimizer = "bobyqa": the default optimizer sometimes fails to converge with
# small datasets. bobyqa (Bound Optimization BY Quadratic Approximation) is more
# robust for small N. optCtrl sets max iterations.
glmer_a <- glmer(
  Is_Positive ~ Country + Technology_Group + (1 | Drug),
  data    = df_glmer_a,
  family  = binomial,             # logistic = binomial family
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(glmer_a))
cat("\n")

# --- Model B: add Severity_Tier (smaller N due to missing values) ------------
#
# Severity_Tier is missing for ~67% of records. Model B filters to the subset
# where it's available, accepting the reduced N. Useful to check whether severity
# still has an effect even after accounting for country and technology.

df_glmer_b <- df_glmer_a |>
  filter(!is.na(Severity_Tier)) |>
  mutate(Severity_Tier = factor(Severity_Tier)) |>
  droplevels()   # removes any unused levels from the filtered subset

# Check how many levels each variable has in the Model B subset.
# If a variable drops to 1 level (e.g. only Norway remains), we can't estimate
# a between-category contrast → exclude it from the formula.
level_check_b <- c(
  Country          = nlevels(df_glmer_b$Country),
  Severity_Tier    = nlevels(df_glmer_b$Severity_Tier),
  Technology_Group = nlevels(df_glmer_b$Technology_Group),
  Evidence_Type    = nlevels(df_glmer_b$Evidence_Type)
)
cat("Model B — N =", nrow(df_glmer_b), "records after Severity_Tier filter\n")
cat("Levels per variable:\n"); print(level_check_b); cat("\n")

# Include only variables with ≥2 levels, AND also remove Evidence_Type explicitly
# (the same separation issue that affects Model A persists in the smaller subset)
valid_b <- names(level_check_b[level_check_b >= 2])
valid_b  <- setdiff(valid_b, "Evidence_Type")    # always exclude Evidence_Type

# Build the formula dynamically from valid_b — prevents errors if a country
# drops out of the filtered subset. paste() + as.formula() constructs the
# formula string programmatically.
formula_b <- as.formula(paste(
  "Is_Positive ~", paste(valid_b, collapse = " + "), "+ (1 | Drug)"
))
cat("Formula used:", deparse(formula_b), "\n\n")

glmer_b <- glmer(
  formula_b,
  data    = df_glmer_b,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(glmer_b))
cat("\n")

# =============================================================================
# SECTION 6: Forest plot — Model A odds ratios
# =============================================================================
#
# ODDS RATIO INTERPRETATION:
# OR = 2.0 → country/technology is associated with 2× higher odds of approval
# OR = 0.5 → associated with 50% lower odds of approval
# OR = 1.0 → no difference from reference
# Reference: Norway (for country), Cell therapy (for technology)
#
# Note: ORs from logistic regression ≠ probability ratios unless events are rare.
# With approval rates ranging from 40–80%, ORs overstate relative risks.
# This is normal for logistic regression and doesn't invalidate the direction
# or significance of associations.
#
# WHY WALD CIs (not likelihood profile)?
# conf.method = "Wald" computes CIs using the quadratic approximation:
#   estimate ± 1.96 × SE
# This is fast but slightly less accurate than profiling, especially for ORs
# far from 1. Profiling is computationally expensive and adds little given
# the wide CIs we already have at N~170.

glmer_a_tidy <- tidy(glmer_a, effects = "fixed", exponentiate = TRUE, conf.int = TRUE,
                     conf.method = "Wald") |>
  filter(term != "(Intercept)") |>    # intercept is not meaningful for forest plots
  mutate(
    # Strip the variable name prefix from factor level labels
    label = term |>
      str_remove("^Country") |>
      str_remove("^Technology_Group") |>
      str_remove("^Evidence_Type") |>
      str_remove("^Severity_Tier"),
    # Assign colour group for the forest plot
    group = case_when(
      str_detect(term, "^Country")          ~ "Country",
      str_detect(term, "^Technology_Group") ~ "Technology",
      str_detect(term, "^Evidence_Type")    ~ "Evidence",
      str_detect(term, "^Severity_Tier")    ~ "Severity",
      TRUE                                  ~ "Other"
    ),
    # Significance stars for quick visual scanning — but ALWAYS check actual p-value
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE             ~ ""
    )
  )

print(glmer_a_tidy |> select(label, group, estimate, conf.low, conf.high, p.value, significance))
cat("\n")

p_forest_glmer <- ggplot(glmer_a_tidy,
       aes(x = estimate, y = reorder(label, estimate), colour = group)) +
  # Reference line at OR = 1 (no effect)
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  # Point estimate + 95% CI whiskers
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.6) +
  # Add significance stars to the right of the CI
  geom_text(aes(label = significance, x = conf.high * 1.05),
            hjust = 0, size = 4, colour = "black") +
  # Log scale: ORs are multiplicative, so 0.5 and 2.0 should be equidistant from 1
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25)) +
  scale_colour_manual(
    values = c(Country = "#1B7EC2", Technology = "#3DAA6D",
               Evidence = "#E8A838", Severity = "#9B59B6")
  ) +
  labs(
    title    = "Mixed-effects logistic regression: odds of positive HTA decision",
    subtitle = "OR > 1 = higher odds of approval; reference: Norway, Cell therapy\nNote: 'Other biologic' and Evidence_Type excluded (quasi-complete separation in small N)\nN=171 records, 50 drugs; drug random intercept; Wald 95% CI",
    x        = "Odds ratio (log scale)",
    y        = NULL,
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(p_forest_glmer)
ggsave(here("figures", "glmer_forest_plot.png"),
       plot = p_forest_glmer, width = 9, height = 6, dpi = 150)

# =============================================================================
# SECTION 7: Quick summary tables
# =============================================================================

cat("Overall approval rate by country:\n")
approval_by_country <- df |>
  filter(!is.na(Is_Positive)) |>
  group_by(Country) |>
  summarise(
    n           = n(),
    positive    = sum(Is_Positive),
    negative    = n() - sum(Is_Positive),
    rate_pct    = round(100 * sum(Is_Positive) / n(), 1)
  )
print(approval_by_country)
cat("\n")

cat("Concordance rate (of drugs evaluated in ≥2 countries):\n")
df_wide |>
  filter(n_countries >= 2) |>
  summarise(
    total                = n(),
    concordant           = sum(is_concordant, na.rm = TRUE),
    discordant           = total - concordant,
    concordance_rate_pct = round(100 * concordant / total, 1)
  ) |>
  print()
cat("\n")

# =============================================================================
# SECTION 8: Excel export
# =============================================================================
#
# Prepare clean export versions of the regression results (rounded, renamed)
# before writing to Excel.

# Model A: tidy output, rounded, rename estimate → OR for clarity
glmer_a_export <- glmer_a_tidy |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  select(label, group, OR = estimate, CI_low = conf.low, CI_high = conf.high,
         p.value, significance)

# Model B: same extraction pipeline
glmer_b_tidy <- tidy(glmer_b, effects = "fixed", exponentiate = TRUE, conf.int = TRUE,
                     conf.method = "Wald") |>
  filter(term != "(Intercept)") |>
  mutate(
    label = term |>
      str_remove("^Country") |> str_remove("^Technology_Group") |>
      str_remove("^Evidence_Type") |> str_remove("^Severity_Tier"),
    group = case_when(
      str_detect(term, "^Country")          ~ "Country",
      str_detect(term, "^Technology_Group") ~ "Technology",
      str_detect(term, "^Evidence_Type")    ~ "Evidence",
      str_detect(term, "^Severity_Tier")    ~ "Severity",
      TRUE                                  ~ "Other"
    ),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE             ~ ""
    )
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  select(label, group, OR = estimate, CI_low = conf.low, CI_high = conf.high,
         p.value, significance)

wb <- createWorkbook()

# Sheet 1: Concordance classification per drug
addWorksheet(wb, "Concordance by Drug")
writeData(wb, "Concordance by Drug", df_wide |>
            select(Drug, Technology_Group, Severity_Tier,
                   pos_Norway, pos_Denmark, pos_Sweden,
                   n_countries, n_positive, concordance))

# Sheet 2: Just the discordant drugs (subset for focused review)
addWorksheet(wb, "Discordant Drugs")
writeData(wb, "Discordant Drugs", discordant_drugs)

# Sheet 3: Approval rate by country (simple proportion table)
addWorksheet(wb, "Approval by Country")
writeData(wb, "Approval by Country", approval_by_country)

# Sheet 4: Mixed-effects Model A coefficients (max N, no Evidence_Type)
addWorksheet(wb, "GLMER Model A (Full N)")
writeData(wb, "GLMER Model A (Full N)", glmer_a_export)

# Sheet 5: Mixed-effects Model B coefficients (with Severity_Tier, smaller N)
addWorksheet(wb, "GLMER Model B (Severity)")
writeData(wb, "GLMER Model B (Severity)", glmer_b_tidy)

saveWorkbook(wb, here("outputs", "03_results.xlsx"), overwrite = TRUE)
cat("Results exported to outputs/03_results.xlsx\n")
cat("Concordance analysis complete.\n")
