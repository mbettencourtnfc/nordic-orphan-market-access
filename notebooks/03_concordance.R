# =============================================================================
# 03_concordance.R
# Nordic Orphan Drug Market Access — Cross-Country Concordance
# =============================================================================
# Goal: When the same drug is evaluated in all three countries (SE, DK, NO),
#       do they agree? What predicts disagreement?
#
# Methods:
#   - Descriptive concordance table (all-agree / partial / all-disagree)
#   - Mixed-effects logistic regression: Is_Positive ~ fixed effects + (1|Drug)
#     Drug as random effect = accounts for the fact that the same drug appears
#     in multiple countries (observations are not independent)
#
# Input:  data/nordic_orphan_clean.csv
# Output: outputs/figures/concordance_heatmap.png
#         outputs/figures/mixed_effects_coefs.png
# =============================================================================

library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(here)

# -----------------------------------------------------------------------------
# 1. Load clean data
# -----------------------------------------------------------------------------
df <- read.csv(here("data", "nordic_orphan_clean.csv"))

# -----------------------------------------------------------------------------
# 2. Concordance summary
# -----------------------------------------------------------------------------
# For drugs evaluated in all 3 countries: classify as
#   "full agreement positive", "full agreement negative", "split decision"
#
# concordance <- df %>%
#   group_by(Drug) %>%
#   filter(n() == 3) %>%   # only drugs with all 3 countries
#   summarise(
#     n_positive = sum(Is_Positive),
#     agreement  = case_when(
#       n_positive == 3 ~ "All positive",
#       n_positive == 0 ~ "All negative",
#       TRUE            ~ "Split"
#     )
#   )
# table(concordance$agreement)

# -----------------------------------------------------------------------------
# 3. Mixed-effects logistic regression
# -----------------------------------------------------------------------------
# Why mixed effects? Each drug appears in up to 3 rows (one per country).
# Treating them as independent would inflate our confidence — drugs that are
# "easy" to approve tend to be approved in all countries. The random effect
# on Drug absorbs this drug-level variation.
#
# model <- glmer(
#   Is_Positive ~ Country + Severity_Tier + Technology_Group + Evidence_Type
#                + (1 | Drug),
#   data   = df,
#   family = binomial
# )
# summary(model)

# -----------------------------------------------------------------------------
# 4. Visualise country agreement patterns
# -----------------------------------------------------------------------------
# Heatmap: drugs on y-axis, countries on x-axis, colour = positive/negative
#
# df_wide <- df %>%
#   select(Drug, Country, Is_Positive) %>%
#   pivot_wider(names_from = Country, values_from = Is_Positive)
#
# TODO: plot heatmap with ggplot2

cat("Concordance analysis skeleton ready — fill in after data cleaning step.\n")
