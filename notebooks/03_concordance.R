# =============================================================================
# 03_concordance.R
# Nordic Orphan Drug Market Access — Cross-Country Concordance Analysis
# =============================================================================
# Research question: When the same orphan drug is evaluated in multiple Nordic
# countries, do they reach the same decision? What predicts agreement/disagreement?
#
# Approach:
#   (A) Descriptive concordance — reshape to wide format, classify each drug as
#       fully concordant, partially concordant, or discordant across countries
#   (B) Mixed-effects logistic regression — models what drives Is_Positive
#       while accounting for the fact that the same drug is assessed in multiple
#       countries (drug as random effect to absorb drug-level confounding)
#
# ⚠ LIMITATIONS:
#   - N=51 drugs; many evaluated in only 1–2 countries → concordance analysis
#     is most meaningful for the subset evaluated in all 3
#   - Mixed-effects model has low power; interpret ORs as exploratory
#   - Severity_Tier missing for ~70% of records → large N drop in models with it
#   - Concordance is for FINAL outcome only; countries may agree in direction but
#     differ in restrictions, conditions, or timing
#
# Input:  data/nordic_orphan_clean.csv
# Output: figures/concordance_heatmap.png
#         figures/concordance_summary.png
#         figures/glmer_forest_plot.png
#         outputs/03_results.xlsx
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(lme4)       # glmer() for mixed-effects logistic regression
library(broom.mixed) # tidy() for glmer objects
library(scales)
library(forcats)
library(stringr)
library(openxlsx)

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("outputs"), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
df <- read_csv(here("data", "nordic_orphan_clean.csv"), show_col_types = FALSE) |>
  mutate(
    Country          = factor(Country, levels = c("Norway", "Denmark", "Sweden")),
    Technology_Group = factor(Technology_Group),
    Evidence_Type    = factor(Evidence_Type),
    Is_Positive      = as.integer(Is_Positive)
  )

cat("Total records loaded:", nrow(df), "\n")
cat("Drugs:", n_distinct(df$Drug), "\n\n")
# -----------------------------------------------------------------------------
# 2. Reshape to wide format — one row per drug
# -----------------------------------------------------------------------------
# For each drug, we want: outcome in each country + key drug-level features
# Drug-level features taken from the first row (same across countries)

drug_features <- df |>
  group_by(Drug) |>
  slice(1) |>
  ungroup() |>
  select(Drug, Technology_Group, Severity_Tier, Therapeutic_Area)

# Some drugs have multiple rows per country (e.g. initial rejection → later approval).
# For concordance we want CURRENT status: if a drug was ever positively reimbursed,
# count it as positive (max over all decisions per Drug × Country).
df_dedup <- df |>
  filter(!is.na(Is_Positive)) |>
  group_by(Drug, Country) |>
  summarise(Is_Positive = max(Is_Positive), .groups = "drop")

cat("Records after deduplication (one per Drug × Country):", nrow(df_dedup), "\n\n")

# Wide format: one column per country
df_wide <- df_dedup |>
  pivot_wider(
    names_from  = Country,
    values_from = Is_Positive,
    names_prefix = "pos_"
  ) |>
  left_join(drug_features, by = "Drug") |>
  mutate(
    n_countries  = rowSums(!is.na(across(starts_with("pos_")))),
    n_positive   = rowSums(across(starts_with("pos_")), na.rm = TRUE),
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
    is_concordant = concordance != "Discordant" & concordance != "Single country only"
  )

cat("Drugs by number of countries evaluated:\n")
print(table(df_wide$n_countries))
cat("\nConcordance overview:\n")
print(table(df_wide$concordance, useNA = "ifany"))
cat("\n")

# Drugs evaluated in all 3 countries
df_all3 <- df_wide |> filter(n_countries == 3)
cat("Drugs evaluated in all 3 countries:", nrow(df_all3), "\n")
cat("  Fully concordant:", sum(df_all3$is_concordant), "\n")
cat("  Discordant:      ", sum(!df_all3$is_concordant), "\n\n")

# Discordant drugs — full details
discordant_drugs <- df_wide |>
  filter(concordance == "Discordant") |>
  arrange(n_countries, Drug) |>
  select(Drug, Technology_Group, Severity_Tier,
         pos_Norway, pos_Denmark, pos_Sweden,
         n_countries, n_positive, concordance)

cat("Discordant drugs:\n")
print(discordant_drugs)
cat("\n")

# -----------------------------------------------------------------------------
# 3. Concordance summary bar chart
# -----------------------------------------------------------------------------
# Shows how many drugs fall into each concordance category, split by n_countries

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
    n_label       = paste0(n_countries, " countr", ifelse(n_countries == 1, "y", "ies"))
  ) |>
  count(n_label, concordance) |>
  ggplot(aes(x = fct_rev(n_label), y = n, fill = concordance)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = concordance_colours, drop = FALSE) +
  coord_flip() +
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
  guides(fill = guide_legend(nrow = 2))

print(p_concordance)
ggsave(here("figures", "concordance_summary.png"),
       plot = p_concordance, width = 9, height = 5, dpi = 150)

# -----------------------------------------------------------------------------
# 4. Outcome heatmap — all drugs × all countries
# -----------------------------------------------------------------------------
# Sort drugs: first by concordance category (fully concordant positive → negative
# → discordant → single country), then alphabetically within each group

drug_order <- df_wide |>
  mutate(sort_key = case_when(
    concordance == "Concordant — all positive"  ~ 1,
    concordance == "Concordant — both positive" ~ 2,
    concordance == "Discordant"                 ~ 3,
    concordance == "Concordant — both negative" ~ 4,
    concordance == "Concordant — all negative"  ~ 5,
    TRUE                                        ~ 6
  )) |>
  arrange(sort_key, Drug) |>
  pull(Drug)

# Long format for heatmap
df_heatmap <- df |>
  select(Drug, Country, Is_Positive) |>
  mutate(
    outcome = case_when(
      Is_Positive == 1 ~ "Positive",
      Is_Positive == 0 ~ "Negative",
      TRUE             ~ NA_character_
    ),
    Drug    = factor(Drug, levels = rev(drug_order)),
    Country = factor(Country, levels = c("Norway", "Denmark", "Sweden"))
  )

p_heatmap <- ggplot(df_heatmap, aes(x = Country, y = Drug, fill = outcome)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(
    values  = c("Positive" = "#3DAA6D", "Negative" = "#cc4444"),
    na.value = "#e8e8e8",
    name    = "HTA decision"
  ) +
  labs(
    title    = "HTA decision outcomes by drug and country",
    subtitle = "Sorted by concordance: all-positive → discordant → all-negative → not evaluated (grey)",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y      = element_text(size = 7),
    legend.position  = "bottom",
    panel.grid       = element_blank()
  )

print(p_heatmap)
ggsave(here("figures", "concordance_heatmap.png"),
       plot = p_heatmap, width = 7, height = 14, dpi = 150)

# -----------------------------------------------------------------------------
# 5. Mixed-effects logistic regression
# -----------------------------------------------------------------------------
# Why mixed-effects?
#   The same drug is assessed in up to 3 countries → observations are correlated.
#   In 02_time_to_access.R we handled this in Cox models with clustered SEs
#   (marginal model approach). Here we use a random-effects approach: add a
#   random intercept per drug, which absorbs unmeasured drug-level factors
#   (clinical benefit, disease severity, price, evidence quality) and lets the
#   fixed effects estimate country/technology/evidence contributions *net of*
#   those drug characteristics.
#
# Fixed effects: Country, Technology_Group, Evidence_Type
# Random effect: (1|Drug)  — drug-specific intercept
#
# Two models:
#   Model A: Country + Technology_Group + Evidence_Type (max N)
#   Model B: adds Severity_Tier (smaller N due to missingness)
#
# Odds ratio (OR) interpretation:
#   OR > 1 = higher odds of positive decision
#   OR < 1 = lower odds
#   Reference: Norway, Cell therapy, RCT evidence

cat("=== Mixed-effects logistic regression ===\n\n")

# --- Model A: maximum N -------------------------------------------------------
# Two predictors cause quasi-complete separation in this small dataset:
#   (a) "Other biologic" — tiny category where all drugs have the same outcome
#   (b) Evidence_Type "Real-world / registry" — all drugs in this evidence type
#       were approved, so the coefficient inflates to ±Inf
# Both cause degenerate Hessians and uninterpretable CIs.
# Fix: exclude "Other biologic" from the model (just 1–2 drugs, noted below);
#      drop Evidence_Type from Model A (run separately as Model A2).
# This is standard practice for small-N logistic models with perfect predictors.

df_glmer <- df |>
  filter(!is.na(Is_Positive)) |>
  mutate(
    Country          = factor(Country, levels = c("Norway", "Denmark", "Sweden")),
    Technology_Group = factor(Technology_Group),
    Evidence_Type    = factor(Evidence_Type),
    Drug             = factor(Drug)
  )

# How many "Other biologic" records?
other_bio_n <- sum(df_glmer$Technology_Group == "Other biologic", na.rm = TRUE)
cat("Note: 'Other biologic' has", other_bio_n,
    "records — excluded from regression (complete separation).\n")
cat("Note: Evidence_Type 'Real-world/registry' also causes separation — dropped from Model A.\n")
cat("      Both are reported in the descriptive tables instead.\n\n")

df_glmer_a <- df_glmer |>
  filter(Technology_Group != "Other biologic") |>
  mutate(Technology_Group = droplevels(Technology_Group))

cat("Model A — Country + Technology_Group + (1|Drug)\n")
cat("N =", nrow(df_glmer_a), "records,", n_distinct(df_glmer_a$Drug), "drugs\n\n")

glmer_a <- glmer(
  Is_Positive ~ Country + Technology_Group + (1 | Drug),
  data    = df_glmer_a,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(glmer_a))
cat("\n")

# --- Model B: add Severity_Tier -----------------------------------------------
df_glmer_b <- df_glmer_a |>
  filter(!is.na(Severity_Tier)) |>
  mutate(Severity_Tier = factor(Severity_Tier)) |>
  droplevels()

# Check levels
level_check_b <- c(
  Country          = nlevels(df_glmer_b$Country),
  Severity_Tier    = nlevels(df_glmer_b$Severity_Tier),
  Technology_Group = nlevels(df_glmer_b$Technology_Group),
  Evidence_Type    = nlevels(df_glmer_b$Evidence_Type)
)
cat("Model B — N =", nrow(df_glmer_b), "records after Severity_Tier filter\n")
cat("Levels per variable:\n"); print(level_check_b); cat("\n")

# Also exclude Evidence_Type from Model B — same quasi-complete separation issue
# in the 57-record subset (all real-world evidence records approved → unestimable)
valid_b <- names(level_check_b[level_check_b >= 2])
valid_b  <- setdiff(valid_b, "Evidence_Type")
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

# -----------------------------------------------------------------------------
# 6. Forest plot — Model A odds ratios
# -----------------------------------------------------------------------------
glmer_a_tidy <- tidy(glmer_a, effects = "fixed", exponentiate = TRUE, conf.int = TRUE,
                     conf.method = "Wald") |>
  filter(term != "(Intercept)") |>
  mutate(
    label = term |>
      str_remove("^Country") |>
      str_remove("^Technology_Group") |>
      str_remove("^Evidence_Type") |>
      str_remove("^Severity_Tier"),
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
  )

print(glmer_a_tidy |> select(label, group, estimate, conf.low, conf.high, p.value, significance))
cat("\n")

p_forest_glmer <- ggplot(glmer_a_tidy,
       aes(x = estimate, y = reorder(label, estimate), colour = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.6) +
  geom_text(aes(label = significance, x = conf.high * 1.05),
            hjust = 0, size = 4, colour = "black") +
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

# -----------------------------------------------------------------------------
# 7. Approval rate by country — quick summary table
# -----------------------------------------------------------------------------
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
    total       = n(),
    concordant  = sum(is_concordant, na.rm = TRUE),
    discordant  = total - concordant,
    concordance_rate_pct = round(100 * concordant / total, 1)
  ) |>
  print()
cat("\n")

# -----------------------------------------------------------------------------
# 8. Excel export
# -----------------------------------------------------------------------------
glmer_a_export <- glmer_a_tidy |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  select(label, group, OR = estimate, CI_low = conf.low, CI_high = conf.high,
         p.value, significance)

glmer_b_tidy <- tidy(glmer_b, effects = "fixed", exponentiate = TRUE, conf.int = TRUE) |>
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

addWorksheet(wb, "Concordance by Drug")
writeData(wb, "Concordance by Drug", df_wide |>
            select(Drug, Technology_Group, Severity_Tier,
                   pos_Norway, pos_Denmark, pos_Sweden,
                   n_countries, n_positive, concordance))

addWorksheet(wb, "Discordant Drugs")
writeData(wb, "Discordant Drugs", discordant_drugs)

addWorksheet(wb, "Approval by Country")
writeData(wb, "Approval by Country", approval_by_country)

addWorksheet(wb, "GLMER Model A (Full N)")
writeData(wb, "GLMER Model A (Full N)", glmer_a_export)

addWorksheet(wb, "GLMER Model B (Severity)")
writeData(wb, "GLMER Model B (Severity)", glmer_b_tidy)

saveWorkbook(wb, here("outputs", "03_results.xlsx"), overwrite = TRUE)
cat("Results exported to outputs/03_results.xlsx\n")
cat("Concordance analysis complete.\n")
