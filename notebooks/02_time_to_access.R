# =============================================================================
# 02_time_to_access.R
# Nordic Orphan Drug Market Access — Time-to-Access Survival Analysis
# =============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Analyses how long it takes for Nordic HTA bodies to make a POSITIVE
# reimbursement decision after EMA marketing authorisation. Uses survival
# analysis — originally developed for clinical time-to-event data, but works
# equally well for policy/market access timelines.
#
# WHY SURVIVAL ANALYSIS?
# ----------------------
# We can't just calculate average time to positive decision because many drugs
# received NEGATIVE decisions. Those drugs might eventually be re-submitted and
# approved — we don't know. In survival terms, they are CENSORED: we know they
# hadn't received access by the time we recorded the data, but we can't say
# they'll never receive access. A simple average would UNDERESTIMATE true time
# to access if we ignored negative decisions, or throw away data if we dropped
# them. Survival analysis uses censored observations without discarding them.
#
# HOW SURVIVAL ANALYSIS IS SET UP HERE
# -------------------------------------
# Time:    time_to_decision_years = years from EMA first MA to HTA decision
# Event:   Is_Positive == TRUE = drug received reimbursement (the "event")
# Censored: Is_Positive == FALSE = drug rejected/not reimbursed by the data
#           collection date; we flag these as event = 0 (censored)
#
# The "survival function" S(t) = probability of still being UNAPPROVED at time t.
# We plot 1 - S(t) = cumulative probability of HAVING been approved by time t.
# This reads more intuitively: "by year 3, ~60% of Swedish drugs were approved."
#
# ASSUMPTIONS
# -----------
# [A1] Negative decisions are treated as right-censored, not as permanent
#      failures. In reality some rejected drugs were resubmitted and approved
#      later. This simplification may slightly underestimate access rates.
# [A2] The same drug appearing in Norway, Denmark, AND Sweden creates 3 rows.
#      These observations are NOT independent — the same drug tends to have
#      similar HTA outcomes everywhere. Cox models handle this via clustered
#      standard errors (see Model A below).
# [A3] The proportional hazards assumption: Cox PH assumes the ratio of hazards
#      (risk of getting approved) between countries stays constant over time.
#      We can't formally test this with N~129, but it's a standard assumption.
# [A4] July-1 date approximation (inherited from 01): adds up to ±6 months noise
#      to time_to_decision_years for ~34 records with year-only dates.
#
# KEY METRICS PRODUCED
# --------------------
# - KM curves: visual picture of how quickly access is achieved by country
# - Log-rank test: formal test of whether countries differ (p-value only)
# - Cox hazard ratios: "how much faster / slower" is each subgroup (with CIs)
# - Forest plot: visual summary of Cox HRs
#
# INTERPRETATION GUIDE
# --------------------
# Hazard ratio (HR) > 1 = faster time to positive decision = BETTER access
# Hazard ratio (HR) < 1 = slower time to positive decision = WORSE access
# Reference category: Norway + Biologics + RCT evidence
# (All other groups compared to these three baselines)
#
# Input:  data/nordic_orphan_clean.csv (from 01_data_cleaning.R)
# Output: figures/km_by_country.png
#         figures/cox_forest_plot.png
#         figures/time_distribution_by_country.png
#         figures/approval_rate_by_technology.png
#         outputs/02_results.xlsx
# =============================================================================

# Load packages ---------------------------------------------------------------
# Each package has a specific role; we don't load tidyverse as a bundle
# because loading only what we need makes dependencies explicit.
library(survival)   # core survival analysis: Surv(), survfit(), coxph(), survdiff()
library(survminer)  # ggsurvplot() — ggplot2-based KM curve visualisation
library(dplyr)      # data manipulation: filter, mutate, group_by, summarise
library(readr)      # read_csv() for the clean dataset
library(ggplot2)    # manual forest plot and distribution plots
library(here)       # here("figures", "...") = project-relative paths
library(broom)      # tidy() turns model output into a clean data frame
library(scales)     # percent_format() for the approval-rate plot
library(stringr)    # str_remove() for cleaning factor level names in plots
library(openxlsx)   # createWorkbook() / saveWorkbook() for Excel export

# Create output folders (safe to run even if they already exist)
# showWarnings = FALSE suppresses the "directory already exists" message
dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)

# =============================================================================
# SECTION 1: Load and prepare data for survival analysis
# =============================================================================
df <- read_csv(here("data", "nordic_orphan_clean.csv"), show_col_types = FALSE)

# Filter to records suitable for survival analysis:
# (a) time_to_decision_years must not be NA (need both EMA and HTA dates)
# (b) time >= 0 only (5 records where HTA date < EMA date are excluded —
#     these are data quality issues, not real negative times)
# (c) Is_Positive must not be NA (we need to know if event occurred)
#
# DECISION: we keep both positive AND negative decisions. Negatives are censored,
# not dropped. This is the whole point of survival analysis — see WHY above.

df_surv <- df |>
  filter(
    !is.na(time_to_decision_years),
    time_to_decision_years >= 0,
    !is.na(Is_Positive)
  ) |>
  mutate(
    # Convert logical TRUE/FALSE to integer 1/0 — required by Surv()
    event          = as.integer(Is_Positive),

    # Set reference level for Country: Norway is reference because it has the
    # most data and is the most commonly cited benchmark in Nordic HTA literature
    Country        = factor(Country, levels = c("Norway", "Denmark", "Sweden")),

    # Factor encoding for Cox covariates
    Technology_Group = factor(Technology_Group),
    Evidence_Type    = factor(Evidence_Type)
  )

# Print data quality summary — these numbers should match what you expect from
# the master data. If something looks off, check 01_data_cleaning.R output.
cat("Records for survival analysis:", nrow(df_surv), "\n")
cat("Events (positive decisions):", sum(df_surv$event), "\n")
cat("Censored (negative):", sum(df_surv$event == 0), "\n\n")

cat("By country:\n")
print(
  df_surv |>
    count(Country, event) |>
    tidyr::pivot_wider(names_from = event, values_from = n,
                       names_prefix = "event_") |>
    rename(positive = event_1, negative_censored = event_0)
)
cat("\n")

# =============================================================================
# SECTION 2: Kaplan-Meier curves by country
# =============================================================================
#
# WHAT IS A KM CURVE?
# A Kaplan-Meier curve estimates the probability of an event NOT having occurred
# yet (survival function) at each time point. Each time an event happens, the
# curve steps down. When a censored observation appears, the person "leaves" the
# risk set but the curve does NOT step down (we just have fewer people left to
# count).
#
# We use fun = "event" to flip it: instead of showing P(not approved yet),
# we show P(approved by time t) = 1 - S(t) = cumulative incidence.
# This means curves going UP and to the RIGHT is GOOD (more approvals faster).

# Create the survival object — this is just a structured pair of (time, event)
# columns that the survival package knows how to handle.
surv_obj <- Surv(time = df_surv$time_to_decision_years, event = df_surv$event)

# Fit one KM curve per country using the formula interface
km_country <- survfit(surv_obj ~ Country, data = df_surv)

cat("Median time to positive HTA decision (years) by country:\n")
print(km_country)
cat("\n")

# Build the KM plot using survminer's ggsurvplot()
# Why ggsurvplot() instead of plain ggplot2?
# Because ggsurvplot handles the "risk table" (how many drugs are still being
# tracked at each time point) automatically — very hard to do in base ggplot2.
#
# pval = TRUE adds the log-rank p-value directly to the plot — quick visual check
# conf.int = TRUE adds 95% confidence bands around each curve
# risk.table = TRUE adds the row of "still at risk" numbers at the bottom

km_plot <- ggsurvplot(
  km_country,
  data        = df_surv,
  fun         = "event",        # flip: show cumulative incidence, not survival
  pval        = TRUE,           # show log-rank p-value on the plot
  conf.int    = TRUE,           # 95% confidence bands
  risk.table  = TRUE,           # risk table below the plot
  risk.table.height = 0.25,     # fraction of total figure height for risk table
  title = paste0(
    "Cumulative probability of positive HTA decision by country\n",
    "(from EMA marketing authorisation; N=", nrow(df_surv), " drug-country pairs)"
  ),
  xlab        = "Years from EMA marketing authorisation",
  ylab        = "Cumulative probability of positive decision",
  legend.title = "Country",
  legend.labs = levels(df_surv$Country),
  palette     = c("#E8A838", "#1B7EC2", "#3DAA6D"), # Norway, Denmark, Sweden
  ggtheme     = theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 11))
)

# print() shows the plot in the VS Code R viewer (requires httpgd)
print(km_plot)

# Save to PNG — overwrites the previous version each time you run this script.
# png() / dev.off() is required for ggsurvplot objects because they contain
# two ggplot panels (curve + risk table) combined with grid.arrange() internally.
# ggsave() only handles single ggplot objects and would fail here.
png(here("figures", "km_by_country.png"), width = 1000, height = 750, res = 120)
print(km_plot)
dev.off()

# =============================================================================
# SECTION 3: Log-rank test for country differences
# =============================================================================
#
# WHAT IS THE LOG-RANK TEST?
# Tests the null hypothesis: "all countries have the same time-to-approval curve."
# It compares observed vs expected events at each event time, summed over time.
# The result is a chi-squared statistic and p-value.
#
# LIMITATION: With N~129 and only 3 groups, we are underpowered (low sample size
# means we might not detect real differences). A p > 0.05 does NOT mean countries
# are the same — it may simply mean we lack evidence with this sample.

logrank_test <- survdiff(surv_obj ~ Country, data = df_surv)

cat("Log-rank test (between-country difference in time to positive decision):\n")
print(logrank_test)
cat("\n")

# Pairwise tests with Bonferroni correction for multiple comparisons.
# Bonferroni: multiplies each p-value by the number of comparisons (3 pairs).
# This is conservative but standard for small N — reduces false positives.
cat("Pairwise log-rank tests (Bonferroni correction):\n")
pairwise_survdiff(Surv(time_to_decision_years, event) ~ Country,
                  data = df_surv, p.adjust.method = "bonferroni") |>
  print()
cat("\n")

# =============================================================================
# SECTION 4: Cox proportional hazards regression
# =============================================================================
#
# WHAT IS COX PH?
# A regression model for time-to-event data. Unlike logistic regression (which
# predicts the PROBABILITY of an event), Cox models predict the RATE of the event
# happening at any given moment (the hazard). The model estimates how each
# predictor shifts the hazard, expressed as a hazard ratio (HR).
#
# HR = 1.5 means: at any given time point, this group's hazard of getting approved
# is 1.5× higher than the reference group = faster / more likely to get approved.
#
# WHY CLUSTERED STANDARD ERRORS?
# The same drug appears in Norway, Denmark, and Sweden → 3 rows in our data.
# These are NOT independent: if Spinraza gets a positive decision in Norway, it's
# more likely to get one in Sweden too (same efficacy, similar disease, similar
# payers). Ordinary Cox assumes independence → underestimates standard errors.
# cluster = Drug tells coxph() to use the sandwich estimator, which produces
# "robust" SEs that account for the within-drug correlation across countries.

# --- Model A: Country + Technology + Evidence (maximum N, excludes Severity_Tier)
# WHY exclude Severity_Tier here?
# Severity_Tier is missing for ~117 / 174 records (~67%). Including it would
# drop 2/3 of our data. Better to run the maximum-N model without it, and
# report Severity_Tier separately in Model B.

cox_a <- coxph(
  Surv(time_to_decision_years, event) ~ Country + Technology_Group + Evidence_Type,
  data    = df_surv,
  cluster = Drug   # robust SEs for within-drug correlation
)

cat("Cox Model A — Country + Technology + Evidence\n")
cat("(clustered SE by drug; N =", nrow(df_surv), ")\n\n")
print(summary(cox_a))
cat("\n")

# --- Model B: add Severity_Tier (accepts that N drops sharply)
# After filtering to records with Severity_Tier, factor levels from OUTSIDE
# that subset may still technically exist as empty levels. For example, if
# "Ultra-severe" only existed in records where Severity_Tier = NA (i.e. none
# of the rows that pass the filter), it would still appear as a level in the
# factor — causing the model to try to estimate a coefficient for a category
# with zero observations → error or NA.
# droplevels() removes any unused factor levels after the filter.

df_surv_b <- df_surv |>
  filter(!is.na(Severity_Tier)) |>
  mutate(Severity_Tier = factor(Severity_Tier, ordered = FALSE)) |>
  droplevels()

cat("Records with Severity_Tier for Model B:", nrow(df_surv_b), "\n")

# Before fitting, check how many distinct levels each variable has in this
# smaller subset. We need ≥2 levels to estimate a coefficient. If a variable
# has only 1 level in the filtered subset (e.g. only Norway remains),
# there's nothing to compare → exclude it.
level_check <- c(
  Country          = nlevels(df_surv_b$Country),
  Severity_Tier    = nlevels(df_surv_b$Severity_Tier),
  Technology_Group = nlevels(df_surv_b$Technology_Group),
  Evidence_Type    = nlevels(df_surv_b$Evidence_Type)
)
cat("Levels per variable in Model B subset:\n")
print(level_check)

valid_vars <- names(level_check[level_check >= 2])
cat("Variables with ≥2 levels (included in Model B):", paste(valid_vars, collapse = ", "), "\n\n")

# Require at least 20 records and 2 valid predictors — below this, Cox estimates
# become unstable and uninterpretable. This check prevents silent bad output.
if (nrow(df_surv_b) >= 20 && length(valid_vars) >= 2) {
  formula_b <- as.formula(
    paste("Surv(time_to_decision_years, event) ~", paste(valid_vars, collapse = " + "))
  )
  cox_b <- coxph(formula_b, data = df_surv_b, cluster = Drug)
  cat("Cox Model B — with Severity_Tier (N =", nrow(df_surv_b), ")\n\n")
  print(summary(cox_b))
  cat("\n")
} else {
  cat("Insufficient data for Model B — skipping.\n\n")
  cox_b <- NULL
}

# =============================================================================
# SECTION 5: Visualisations
# =============================================================================

# --- 5a. Forest plot: Cox Model A hazard ratios ----------------------------
#
# WHAT IS A FOREST PLOT?
# A forest plot shows the point estimates and confidence intervals for each
# predictor in a regression model. The vertical dashed line at HR = 1 is the
# "no effect" reference. Points to the RIGHT of 1 mean faster access; to the
# LEFT mean slower access.
#
# WHY NOT USE ggforest() FROM survminer?
# ggforest() is a convenience function that plots Cox model results, but it does
# NOT support the cluster = argument. Passing our clustered Cox model causes an
# error. We replicate it manually using broom::tidy() to extract results and
# ggplot2 to plot them.
#
# NOTE ON CONFIDENCE INTERVALS:
# tidy() produces NAIVE (unclustered) CIs. The point estimates (HRs) are the
# same as the clustered model, but the CIs are slightly too narrow because they
# don't account for the within-drug correlation. This is a known limitation of
# tidy() with clustered Cox models — exact clustered CIs require manual extraction.
# This is acceptable for exploratory / portfolio purposes; note it in interpretation.

forest_data <- tidy(cox_a, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    # Strip the variable prefix R adds to factor level names.
    # e.g. "CountryDenmark" → "Denmark", "Technology_GroupGene therapy" → "Gene therapy"
    label = term |>
      str_remove("^Country")           |>
      str_remove("^Technology_Group")  |>
      str_remove("^Evidence_Type")     ,
    # Colour-code each row by the variable it belongs to
    group = case_when(
      str_detect(term, "^Country")          ~ "Country",
      str_detect(term, "^Technology_Group") ~ "Technology",
      str_detect(term, "^Evidence_Type")    ~ "Evidence"
    )
  )

forest_a <- ggplot(forest_data,
       aes(x = estimate, y = reorder(label, estimate), colour = group)) +
  # Vertical reference line at HR = 1 (no effect)
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  # Point for HR estimate + whiskers for 95% CI
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.6) +
  # Log scale because HRs are multiplicative (0.5 and 2.0 are equally "far" from 1)
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  scale_colour_manual(
    values = c(Country = "#1B7EC2", Technology = "#3DAA6D", Evidence = "#E8A838")
  ) +
  labs(
    title    = "Cox Model A: hazard ratios for time to positive HTA decision",
    subtitle = "HR > 1 = faster access; reference: Norway, Biologics, RCT evidence\nN=129 drug-country pairs, clustered SE by drug",
    x        = "Hazard ratio (log scale)",
    y        = NULL,
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(forest_a)
# ggsave() works here because forest_a is a single standard ggplot2 object
ggsave(here("figures", "cox_forest_plot.png"), plot = forest_a, width = 9, height = 6, dpi = 150)

# --- 5b. Time-to-decision distribution by country --------------------------
#
# Violin + boxplot + jitter: three layers of information:
# - Violin: shape of the distribution (is it bimodal? skewed?)
# - Boxplot: median, IQR, outliers
# - Jitter: individual data points, coloured by outcome (positive = black, negative = red)
#
# WHY NOT JUST A BOXPLOT?
# With N~40 per country, individual points matter. A boxplot hides them.
# Showing points coloured by outcome reveals, for example, whether the long-time
# records are mostly negative decisions (rejected late) or approvals (slow process).

p_time <- ggplot(df_surv, aes(x = Country, y = time_to_decision_years, fill = Country)) +
  geom_violin(alpha = 0.35, colour = NA) +       # shape of distribution
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +  # summary stats
  geom_jitter(aes(colour = factor(event)), width = 0.12, size = 1.8, alpha = 0.8) + # points
  scale_fill_manual(values  = c(Norway = "#E8A838", Denmark = "#1B7EC2", Sweden = "#3DAA6D"),
                    guide   = "none") +           # suppress fill legend (redundant with x axis)
  scale_colour_manual(values = c("0" = "#cc3333", "1" = "#333333"),
                      labels = c("0" = "Negative / censored", "1" = "Positive decision"),
                      name   = "HTA outcome") +
  labs(
    title    = "Time from EMA authorisation to HTA decision by country",
    subtitle = paste0("N=", nrow(df_surv), " drug-country pairs with valid time data"),
    x        = NULL,
    y        = "Years from EMA marketing authorisation"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_time)
ggsave(here("figures", "time_distribution_by_country.png"), plot = p_time, width = 8, height = 6, dpi = 150)

# --- 5c. Approval rate by country and technology group ---------------------
#
# This is a simpler descriptive plot — not survival analysis, just the raw
# proportion of drugs that received a positive decision, broken down by
# technology type and country.
#
# DECISION: suppress cells with n < 2 to avoid misleading 0% or 100% rates
# from a single drug. A single approved gene therapy doesn't mean "100% of
# gene therapies are approved in Denmark" — it's just one data point.

p_rate <- df |>
  filter(!is.na(Is_Positive), !is.na(Technology_Group)) |>
  mutate(Country = factor(Country, levels = c("Norway", "Denmark", "Sweden"))) |>
  group_by(Country, Technology_Group) |>
  summarise(
    n         = n(),
    approved  = sum(Is_Positive, na.rm = TRUE),
    rate      = approved / n,
    .groups   = "drop"
  ) |>
  filter(n >= 2) |>   # suppress small cells
  ggplot(aes(x = reorder(Technology_Group, rate), y = rate, fill = Country)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  # Label each bar with "approved / total" so readers can see raw numbers,
  # not just percentages which can be misleading at small N
  geom_text(aes(label = paste0(approved, "/", n)),
            position = position_dodge(0.8), hjust = -0.1, size = 2.8) +
  scale_fill_manual(values = c(Norway = "#E8A838", Denmark = "#1B7EC2", Sweden = "#3DAA6D")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.15)) +
  coord_flip() +
  labs(
    title    = "Positive HTA decision rate by technology group and country",
    subtitle = "Cells with <2 observations suppressed",
    x        = NULL,
    y        = "Approval rate"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(p_rate)
ggsave(here("figures", "approval_rate_by_technology.png"), plot = p_rate, width = 9, height = 6, dpi = 150)

# =============================================================================
# SECTION 6: Descriptive tables (computed here, exported below)
# =============================================================================
#
# These tables appear in the Excel output. Computed once and stored as objects
# so we can reuse them without recalculating.

# Time-to-decision summary by country
# q25/q75 = first and third quartile = the range covering the middle 50% of drugs
tbl_country <- df_surv |>
  group_by(Country) |>
  summarise(
    n                  = n(),
    positive_decisions = sum(event),
    negative_censored  = n() - sum(event),
    approval_rate_pct  = round(100 * sum(event) / n(), 1),
    median_years       = round(median(time_to_decision_years), 2),
    q25_years          = round(quantile(time_to_decision_years, 0.25), 2),
    q75_years          = round(quantile(time_to_decision_years, 0.75), 2),
    max_years          = round(max(time_to_decision_years), 2)
  )

# Same summary by technology group — tells us whether certain technology types
# (e.g. gene therapy) tend to take longer or have lower approval rates
tbl_technology <- df_surv |>
  group_by(Technology_Group) |>
  summarise(
    n                  = n(),
    positive_decisions = sum(event),
    approval_rate_pct  = round(100 * sum(event) / n(), 1),
    median_years       = round(median(time_to_decision_years), 2),
    q25_years          = round(quantile(time_to_decision_years, 0.25), 2),
    q75_years          = round(quantile(time_to_decision_years, 0.75), 2)
  ) |>
  arrange(median_years)

cat("Descriptive time-to-decision by country (years):\n")
print(tbl_country)
cat("\n")
cat("Descriptive time-to-decision by technology group (years):\n")
print(tbl_technology)
cat("\n")

# KM estimates table: median survival time per country from the KM fit
# as.data.frame() + rownames_to_column() converts the matrix summary to a tibble
km_summary <- summary(km_country)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("Group") |>
  mutate(Group = str_remove(Group, "Country="))

# Clean Cox Model A results for export
# exponentiate = TRUE converts log-HR to HR (easier to interpret)
# The significance stars are for quick visual scanning — always look at the actual
# p-value and CI, not just the stars, especially with N < 200
cox_a_tidy <- tidy(cox_a, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    term = term |>
      str_remove("^Country") |>
      str_remove("^Technology_Group") |>
      str_remove("^Evidence_Type"),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE             ~ ""
    )
  ) |>
  rename(HR = estimate, CI_low = conf.low, CI_high = conf.high) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# Cox Model B results (conditional on whether Model B ran successfully)
cox_b_tidy <- if (!is.null(cox_b)) {
  tidy(cox_b, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(
      term = term |>
        str_remove("^Country") |>
        str_remove("^Technology_Group") |>
        str_remove("^Evidence_Type"),
      significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE             ~ ""
      )
    ) |>
    rename(HR = estimate, CI_low = conf.low, CI_high = conf.high) |>
    mutate(across(where(is.numeric), \(x) round(x, 4)))
} else {
  # Placeholder row so the Excel sheet isn't empty
  data.frame(note = "Model B not fitted — insufficient data after Severity_Tier filter")
}

# Summarise the log-rank test as a one-row table for the Excel export
# pchisq() converts the chi-squared statistic to a p-value
logrank_p  <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
tbl_logrank <- data.frame(
  test          = "Log-rank (overall)",
  chi_squared   = round(logrank_test$chisq, 4),
  df            = length(logrank_test$n) - 1,
  p_value       = round(logrank_p, 4),
  interpretation = ifelse(logrank_p < 0.05,
                          "Significant difference between countries",
                          "No significant difference (possibly underpowered)")
)

# Approval rate matrix: Country × Technology_Group
# Shows raw fraction, not survival-adjusted. Useful cross-check for the Cox results.
tbl_approval_matrix <- df |>
  filter(!is.na(Is_Positive), !is.na(Technology_Group)) |>
  group_by(Country, Technology_Group) |>
  summarise(n = n(), approved = sum(Is_Positive, na.rm = TRUE),
            rate_pct = round(100 * sum(Is_Positive, na.rm = TRUE) / n(), 1),
            .groups = "drop") |>
  mutate(suppressed = n < 2)   # flag cells below the display threshold

# =============================================================================
# SECTION 7: Excel export — all numerical results to one file
# =============================================================================
#
# WHY EXCEL AND NOT JUST CSV?
# Multiple sheets let us keep all related tables in one file — easier to share
# and read than a folder of separate CSV files. openxlsx doesn't require Java
# (unlike xlsx package), so it's more reliable across operating systems.

wb <- createWorkbook()

# Sheet 1: High-level country comparison
addWorksheet(wb, "Country Summary")
writeData(wb, "Country Summary", tbl_country)

# Sheet 2: Technology group breakdown
addWorksheet(wb, "Technology Summary")
writeData(wb, "Technology Summary", tbl_technology)

# Sheet 3: Kaplan-Meier formal estimates (including median survival, n.events)
addWorksheet(wb, "KM Estimates")
writeData(wb, "KM Estimates", km_summary)

# Sheet 4: Log-rank test result
addWorksheet(wb, "Log-Rank Test")
writeData(wb, "Log-Rank Test", tbl_logrank)

# Sheet 5: Cox Model A coefficients (HR, CI, p)
addWorksheet(wb, "Cox Model A")
writeData(wb, "Cox Model A", cox_a_tidy)

# Sheet 6: Cox Model B coefficients (with Severity_Tier, smaller N)
addWorksheet(wb, "Cox Model B")
writeData(wb, "Cox Model B", cox_b_tidy)

# Sheet 7: Approval rates cross-tabulated by Country × Technology
addWorksheet(wb, "Approval by Tech x Country")
writeData(wb, "Approval by Tech x Country", tbl_approval_matrix)

# overwrite = TRUE replaces the file each run (no versioned accumulation)
saveWorkbook(wb, here("outputs", "02_results.xlsx"), overwrite = TRUE)
cat("Results exported to outputs/02_results.xlsx\n")
cat("Time-to-access analysis complete.\n")
