# =============================================================================
# 02_time_to_access.R
# Nordic Orphan Drug Market Access — Time-to-Access Analysis
# =============================================================================
# Research question: After EMA marketing authorisation, how long do Nordic
# HTA bodies take to make a positive reimbursement decision? Do countries differ?
#
# Survival framework:
#   Time:  time_to_decision_years (EMA MA date → HTA decision date)
#   Event: Is_Positive == TRUE (drug received positive/partial reimbursement)
#   Censored: negative decisions (drug could be resubmitted; outcome not yet
#             final from access perspective)
#
# ⚠ LIMITATIONS:
#   - N=51 drugs × 3 countries = 174 total records; ~129 with valid time data
#   - Small N means wide CIs and low power — interpret descriptively
#   - 34 records use July-1 date approximation (Decision_Year only)
#   - Negative decisions treated as censored (simplification)
#   - Same drug across 3 countries creates correlated observations →
#     Cox models use clustered SEs to account for within-drug correlation
#
# Input:  data/nordic_orphan_clean.csv
# Output: figures/km_by_country.png
#         figures/cox_forest_plot.png
# =============================================================================

library(survival)
library(survminer)
library(dplyr)
library(readr)
library(ggplot2)
library(here)
library(broom)      # tidy() extracts model coefficients cleanly
library(scales)     # percent_format() for approval rate chart
library(stringr)    # str_remove() for cleaning label names
library(openxlsx)   # Excel export

# Figures folder — created once, plots always overwrite previous versions
dir.create(here("figures"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load and prepare
# -----------------------------------------------------------------------------
df <- read_csv(here("data", "nordic_orphan_clean.csv"), show_col_types = FALSE)

df_surv <- df |>
  filter(
    !is.na(time_to_decision_years),
    time_to_decision_years >= 0,  # exclude 5 records where HTA pre-dates EMA MA
    !is.na(Is_Positive)
  ) |>
  mutate(
    event          = as.integer(Is_Positive),
    Country        = factor(Country, levels = c("Norway", "Denmark", "Sweden")),
    Technology_Group = factor(Technology_Group),
    Evidence_Type  = factor(Evidence_Type)
  )

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

# -----------------------------------------------------------------------------
# 2. Kaplan-Meier curves by country
# -----------------------------------------------------------------------------
surv_obj   <- Surv(time = df_surv$time_to_decision_years, event = df_surv$event)
km_country <- survfit(surv_obj ~ Country, data = df_surv)

cat("Median time to positive HTA decision (years) by country:\n")
print(km_country)
cat("\n")

# Plot cumulative incidence (1 - KM): probability of HAVING been approved by time t
# This reads more naturally as "access probability over time"
km_plot <- ggsurvplot(
  km_country,
  data        = df_surv,
  fun         = "event",      # cumulative incidence = 1 - survival
  pval        = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  risk.table.height = 0.25,
  title = paste0(
    "Cumulative probability of positive HTA decision by country\n",
    "(from EMA marketing authorisation; N=", nrow(df_surv), " drug-country pairs)"
  ),
  xlab        = "Years from EMA marketing authorisation",
  ylab        = "Cumulative probability of positive decision",
  legend.title = "Country",
  legend.labs = levels(df_surv$Country),
  palette     = c("#E8A838", "#1B7EC2", "#3DAA6D"),
  ggtheme     = theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 11))
)

print(km_plot)
# Save to figures (overwrites previous version)
png(here("figures", "km_by_country.png"), width = 1000, height = 750, res = 120)
print(km_plot)
dev.off()

# -----------------------------------------------------------------------------
# 3. Log-rank test
# -----------------------------------------------------------------------------
logrank_test <- survdiff(surv_obj ~ Country, data = df_surv)

cat("Log-rank test (between-country difference in time to positive decision):\n")
print(logrank_test)
cat("\n")

# Pairwise log-rank (Bonferroni corrected)
cat("Pairwise log-rank tests (Bonferroni correction):\n")
pairwise_survdiff(Surv(time_to_decision_years, event) ~ Country,
                  data = df_surv, p.adjust.method = "bonferroni") |>
  print()
cat("\n")

# -----------------------------------------------------------------------------
# 4. Cox proportional hazards
# -----------------------------------------------------------------------------
# HR > 1 = faster time to positive decision (higher "hazard" of being approved)
# HR < 1 = slower / less likely to be approved at any given time point
#
# Model A: country + technology + evidence (maximises N by excluding Severity_Tier,
#          which has 117/174 missing values)
#
# Clustered SE: same drug appears in up to 3 countries → observations correlated

cox_a <- coxph(
  Surv(time_to_decision_years, event) ~ Country + Technology_Group + Evidence_Type,
  data    = df_surv,
  cluster = Drug
)

cat("Cox Model A — Country + Technology + Evidence\n")
cat("(clustered SE by drug; N =", nrow(df_surv), ")\n\n")
print(summary(cox_a))
cat("\n")

# Model B: add Severity_Tier (accepted that N drops sharply)
# After filtering to records with Severity_Tier, some factor levels may
# disappear entirely. droplevels() removes them; then we check which variables
# still have ≥2 levels before including them in the formula.
df_surv_b <- df_surv |>
  filter(!is.na(Severity_Tier)) |>
  mutate(Severity_Tier = factor(Severity_Tier, ordered = FALSE)) |>
  droplevels()

cat("Records with Severity_Tier for Model B:", nrow(df_surv_b), "\n")

# Check which variables retain ≥2 levels in this subset
level_check <- c(
  Country          = nlevels(df_surv_b$Country),
  Severity_Tier    = nlevels(df_surv_b$Severity_Tier),
  Technology_Group = nlevels(df_surv_b$Technology_Group),
  Evidence_Type    = nlevels(df_surv_b$Evidence_Type)
)
cat("Levels per variable in Model B subset:\n")
print(level_check)

# Only include variables with ≥2 levels
valid_vars <- names(level_check[level_check >= 2])
cat("Variables with ≥2 levels (included in Model B):", paste(valid_vars, collapse = ", "), "\n\n")

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

# -----------------------------------------------------------------------------
# 5. Visualisations
# -----------------------------------------------------------------------------

# --- 5a. Forest plot: Cox Model A hazard ratios ----------------------------
# Using broom::tidy() + ggplot2 instead of ggforest() — more reliable with
# clustered models and gives us full control over labels and appearance.
# HRs come from the clustered model (cox_a); CIs are approximate (naive SE)
# since tidy() doesn't expose clustered CIs directly. Note this in interpretation.

forest_data <- tidy(cox_a, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    # Clean up the R-generated factor label names
    label = term |>
      str_remove("^Country")           |>
      str_remove("^Technology_Group")  |>
      str_remove("^Evidence_Type")     ,
    # Group each row for colour coding
    group = case_when(
      str_detect(term, "^Country")          ~ "Country",
      str_detect(term, "^Technology_Group") ~ "Technology",
      str_detect(term, "^Evidence_Type")    ~ "Evidence"
    )
  )

forest_a <- ggplot(forest_data,
       aes(x = estimate, y = reorder(label, estimate), colour = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.6) +
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
ggsave(here("figures", "cox_forest_plot.png"), plot = forest_a, width = 9, height = 6, dpi = 150)

# --- 5b. Time-to-decision distribution by country --------------------------
# Violin + boxplot + jitter: shows spread, not just the median.
# Only records with valid time AND a definitive outcome (positive or negative)
p_time <- ggplot(df_surv, aes(x = Country, y = time_to_decision_years, fill = Country)) +
  geom_violin(alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = factor(event)), width = 0.12, size = 1.8, alpha = 0.8) +
  scale_fill_manual(values  = c(Norway = "#E8A838", Denmark = "#1B7EC2", Sweden = "#3DAA6D"),
                    guide   = "none") +
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
  filter(n >= 2) |>  # suppress cells with <2 observations
  ggplot(aes(x = reorder(Technology_Group, rate), y = rate, fill = Country)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
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

# -----------------------------------------------------------------------------
# 7. Descriptive tables (saved to variables for Excel export below)
# -----------------------------------------------------------------------------
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

# KM estimates table
km_summary <- summary(km_country)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("Group") |>
  mutate(Group = str_remove(Group, "Country="))

# Cox Model A tidy results
cox_a_tidy <- tidy(cox_a, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    term        = term |> str_remove("^Country") |>
                    str_remove("^Technology_Group") |> str_remove("^Evidence_Type"),
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

# Cox Model B tidy results (if it ran)
cox_b_tidy <- if (!is.null(cox_b)) {
  tidy(cox_b, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(
      term = term |> str_remove("^Country") |>
               str_remove("^Technology_Group") |> str_remove("^Evidence_Type"),
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
  data.frame(note = "Model B not fitted — insufficient data after Severity_Tier filter")
}

# Log-rank test summary
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

# Approval rates by country × technology
tbl_approval_matrix <- df |>
  filter(!is.na(Is_Positive), !is.na(Technology_Group)) |>
  group_by(Country, Technology_Group) |>
  summarise(n = n(), approved = sum(Is_Positive, na.rm = TRUE),
            rate_pct = round(100 * sum(Is_Positive, na.rm = TRUE) / n(), 1),
            .groups = "drop") |>
  mutate(suppressed = n < 2)

# -----------------------------------------------------------------------------
# 8. Excel export — outputs/02_results.xlsx
# -----------------------------------------------------------------------------
dir.create(here("outputs"), showWarnings = FALSE)

wb <- createWorkbook()

addWorksheet(wb, "Country Summary")
writeData(wb, "Country Summary", tbl_country)

addWorksheet(wb, "Technology Summary")
writeData(wb, "Technology Summary", tbl_technology)

addWorksheet(wb, "KM Estimates")
writeData(wb, "KM Estimates", km_summary)

addWorksheet(wb, "Log-Rank Test")
writeData(wb, "Log-Rank Test", tbl_logrank)

addWorksheet(wb, "Cox Model A")
writeData(wb, "Cox Model A", cox_a_tidy)

addWorksheet(wb, "Cox Model B")
writeData(wb, "Cox Model B", cox_b_tidy)

addWorksheet(wb, "Approval by Tech x Country")
writeData(wb, "Approval by Tech x Country", tbl_approval_matrix)

saveWorkbook(wb, here("outputs", "02_results.xlsx"), overwrite = TRUE)
cat("Results exported to outputs/02_results.xlsx\n")

cat("Time-to-access analysis complete.\n")
