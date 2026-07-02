# =============================================================================
# 00_fetch_ema_ma_dates.R
# Get the EMA first marketing authorisation date for each of our 51 drugs
# =============================================================================
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# To compute "time from EMA approval to HTA decision" we need to know WHEN
# each drug was first authorised by the European Medicines Agency (EMA).
#
# You might think: just download the EMA orphan designation list (COMP file).
# We did — and it didn't work. The "European Commission decision date" in that
# file is the date of the MOST RECENT EC decision, which is often a 2023/2024
# renewal or amendment, not the original 2015/2016 first authorisation. Using
# it gave us negative time values (HTA decisions before EMA approval) for old
# drugs, which is impossible.
#
# The correct source is the EMA medicines data table, which has a
# "Marketing authorisation date" column = the actual date EMA first authorised
# the medicine for sale. This is updated nightly by EMA and freely downloadable.
#
# RUN THIS ONCE before running 01_data_cleaning.R.
# If the EMA file hasn't changed, you don't need to re-run it.
#
# ASSUMPTIONS
# -----------
# [A1] We trust the "Marketing authorisation date" in the EMA medicines table
#      as the ground truth for first EU approval. For drugs with conditional MA
#      or that were later withdrawn, this is still the original grant date.
# [A2] For combination drugs (e.g. "lumacaftor/ivacaftor"), the EMA file lists
#      the full combination string. Our master data may spell it differently.
#      We handle this with fuzzy matching and manual overrides.
# [A3] Two drugs in our dataset have NO EMA central authorisation:
#      - delandistrogene moxeparvovec (Elevidys): FDA approved 2023, EMA review
#        still ongoing as of mid-2024. Left as NA.
#      - birch triterpenes: Norwegian-specific product, no EMA central MA. NA.
#      These two are excluded from any time-to-decision analysis automatically
#      because their EMA dates are NA.
#
# Output: data/raw/ema_ma_dates.csv
#         Columns: inn_clean, ema_ma_date, orphan_flag
# =============================================================================

library(readxl)    # read EMA Excel file
library(dplyr)     # data manipulation
library(stringr)   # string cleaning and matching
library(here)      # relative file paths (always points to project root)
library(lubridate) # date parsing

# -----------------------------------------------------------------------------
# 1. Download the EMA medicines Excel
# -----------------------------------------------------------------------------
# EMA publishes an Excel file of all authorised medicines at this URL.
# It is updated nightly. The file is ~2 MB.
# We save it locally so we can work with it without re-downloading each time.

ema_url   <- "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx"
dest_file <- here("data", "raw", "ema_medicines_raw.xlsx")

cat("Downloading EMA medicines data from ema.europa.eu...\n")
cat("(~2 MB — takes a few seconds)\n")
download.file(ema_url, destfile = dest_file, mode = "wb", quiet = FALSE)
cat("Download complete.\n\n")

# -----------------------------------------------------------------------------
# 2. Read and parse the EMA file
# -----------------------------------------------------------------------------
# IMPORTANT: The real column headers are on ROW 9, not row 1.
# Rows 1–8 are metadata/title rows. skip = 8 tells R to skip those.
# If you open the file in Excel and see the headers on row 9, that's why.

ema_raw <- read_excel(dest_file, skip = 8)

cat("EMA file columns (first 10):\n")
print(names(ema_raw)[1:10])
cat("\n")

# Extract only the columns we need, for human medicines that are currently authorised
# DECISION: We filter to Category == "Human" and Medicine status == "Authorised"
# This excludes veterinary medicines and withdrawn/refused medicines.
# Drugs with conditional MA or exceptional circumstances ARE included —
# their authorisation date is still their first authorisation date.

ema_ma_dates <- ema_raw |>
  filter(Category == "Human", `Medicine status` == "Authorised") |>
  select(
    inn_ema     = `International non-proprietary name (INN) / common name`,
    ma_date_raw = `Marketing authorisation date`,
    orphan_flag = `Orphan medicine`
  ) |>
  filter(!is.na(inn_ema), !is.na(ma_date_raw)) |>
  mutate(
    # Standardise INN to lowercase with no leading/trailing whitespace.
    # This is our join key — both sides need to look the same for matching to work.
    inn_clean = str_to_lower(str_trim(inn_ema)),

    # The date column comes out as numeric (Excel serial date format) OR
    # as a character string depending on your R/readxl version.
    # We handle both cases here.
    # Excel stores dates as days since 1899-12-30, so we convert with that origin.
    ema_ma_date = case_when(
      is.numeric(ma_date_raw) ~ as_date(as.numeric(ma_date_raw), origin = "1899-12-30"),
      TRUE                    ~ dmy(as.character(ma_date_raw))
    )
  ) |>
  filter(!is.na(ema_ma_date)) |>

  # ASSUMPTION [A1a]: If a drug has multiple authorisation records (e.g. original
  # MA + line extension), we take the EARLIEST date — the true first authorisation.
  group_by(inn_clean) |>
  slice_min(ema_ma_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(inn_clean, ema_ma_date, orphan_flag)

cat("EMA MA dates extracted:", nrow(ema_ma_dates), "unique authorised human medicines\n")
cat("Date range:", format(min(ema_ma_dates$ema_ma_date)), "to",
    format(max(ema_ma_dates$ema_ma_date)), "\n\n")

# -----------------------------------------------------------------------------
# 3. Check coverage against our 51-drug dataset
# -----------------------------------------------------------------------------
master   <- read_excel(here("data", "raw", "nordic_orphan_master.xlsx"), sheet = "Master Data")
our_inns <- str_to_lower(str_trim(unique(master$INN)))

matched   <- our_inns[our_inns %in% ema_ma_dates$inn_clean]
unmatched <- our_inns[!our_inns %in% ema_ma_dates$inn_clean]

cat("Direct INN matches:", length(matched), "/", length(our_inns), "\n")
if (length(unmatched) > 0) {
  cat("Unmatched INNs (will attempt fuzzy matching):\n")
  print(unmatched)
}
cat("\n")

# -----------------------------------------------------------------------------
# 4. Fuzzy matching for unmatched drugs
# -----------------------------------------------------------------------------
# Some drugs don't match because:
#   (a) Combination drugs: our dataset says "lumacaftor/ivacaftor" but EMA
#       lists each component separately, OR stores the combo differently.
#   (b) Gene therapy brand names: EMA sometimes uses a slightly different INN
#       spelling or the drug is listed under its first word only.
#
# Strategy:
#   Step (a): For combination drugs (containing "/"), try matching on just
#             the FIRST component. E.g. "lumacaftor/ivacaftor" → try "lumacaftor".
#             If that matches an EMA entry, use that date.
#             CAVEAT: this works only if the combination drug was approved on the
#             same date as the first component — which isn't always true. For the
#             drugs in our dataset this was confirmed by manual check.
#
#   Step (b): For non-combination drugs, try matching the FIRST WORD of the INN.
#             E.g. "betibeglogene autotemcel" → try "betibeglogene".
#             This catches truncation/spacing differences.

# (a) Combination drugs
combo_unmatched <- unmatched[grepl("/", unmatched)]

combo_fix <- data.frame(inn_master = combo_unmatched) |>
  mutate(
    first_component = str_to_lower(str_trim(str_extract(inn_master, "^[^/]+")))
  ) |>
  left_join(
    ema_ma_dates |> select(first_component = inn_clean, ema_ma_date, orphan_flag),
    by = "first_component"
  ) |>
  filter(!is.na(ema_ma_date)) |>
  select(inn_clean = inn_master, ema_ma_date, orphan_flag)

cat("Combination drug matches (first-component strategy):", nrow(combo_fix), "\n")
if (nrow(combo_fix) > 0) print(combo_fix)
cat("\n")

# (b) First-word matching for remaining single-component unmatched drugs
still_unmatched <- unmatched[!unmatched %in% combo_fix$inn_clean]
still_unmatched <- still_unmatched[!grepl("/", still_unmatched)]

fuzzy_fix <- lapply(still_unmatched, function(our_inn) {
  first_word <- str_extract(our_inn, "^\\S+")
  hit <- ema_ma_dates |>
    filter(str_detect(inn_clean, fixed(first_word, ignore_case = TRUE))) |>
    slice_min(ema_ma_date, n = 1, with_ties = FALSE)
  if (nrow(hit) == 0) return(NULL)
  tibble(inn_clean = our_inn, ema_ma_date = hit$ema_ma_date,
         orphan_flag = hit$orphan_flag, matched_on = hit$inn_clean)
}) |> bind_rows()

cat("First-word fuzzy matches:", nrow(fuzzy_fix), "\n")
if (nrow(fuzzy_fix) > 0) print(fuzzy_fix)
cat("\n")

# -----------------------------------------------------------------------------
# 5. Manual overrides for drugs not caught by any automatic matching
# -----------------------------------------------------------------------------
# These drugs required manual lookup on individual EMA medicine pages.
# Source for each: ema.europa.eu/en/medicines/human/EPAR/<medicine-name>
#
# Why they needed manual handling:
# - Combination drugs (Orkambi, Kaftrio): EMA stores them as the full combination
#   string with spaces around slashes, which our fuzzy match missed.
# - Gene therapies with conditional MA (Zynteglo/betibeglogene, Libmeldy, Skysona):
#   These show as "Withdrawn" in the EMA file because the MA was later withdrawn
#   or replaced, so they were filtered out by `Medicine status == "Authorised"`.
#   We still want their ORIGINAL authorisation date for computing time-to-HTA.
# - Fidanacogene elaparvovec (Beqvez): very recent approval (Sept 2024), INN
#   spelling mismatch.
# - Ataluren (Translarna): conditional MA granted 2014, status changed 2023.
#
# Two drugs are left as NA intentionally (see ASSUMPTION [A3] in header):
#   - delandistrogene moxeparvovec (Elevidys)
#   - birch triterpenes

manual_overrides <- tribble(
  ~inn_clean,                           ~ema_ma_date,           ~orphan_flag,
  "lumacaftor/ivacaftor",               as_date("2015-11-19"),  "Yes",   # Orkambi
  "elexacaftor/tezacaftor/ivacaftor",   as_date("2020-08-21"),  "Yes",   # Kaftrio
  "betibeglogene autotemcel",           as_date("2019-08-29"),  "Yes",   # Zynteglo (conditional, later withdrawn)
  "atidarsagene autotemcel",            as_date("2020-12-16"),  "Yes",   # Libmeldy
  "elivaldogene autotemcel",            as_date("2021-07-26"),  "Yes",   # Skysona
  "fidanacogene elaparvovec",           as_date("2024-09-19"),  "Yes",   # Beqvez
  "ataluren",                           as_date("2014-07-31"),  "Yes",   # Translarna (conditional MA)
) |> select(inn_clean, ema_ma_date, orphan_flag)

cat("Manual overrides:", nrow(manual_overrides), "drugs\n\n")

# -----------------------------------------------------------------------------
# 6. Combine all sources and save
# -----------------------------------------------------------------------------
# Priority order for deduplication: EMA direct match > combo fix > fuzzy > manual
# In practice most drugs only appear in one source.
# distinct(.keep_all = TRUE) keeps the FIRST occurrence — so direct EMA matches
# win over manual overrides if both exist.

ema_ma_dates_full <- bind_rows(
  ema_ma_dates |> select(inn_clean, ema_ma_date, orphan_flag),
  combo_fix,
  fuzzy_fix |> select(inn_clean, ema_ma_date, orphan_flag),
  manual_overrides
) |>
  distinct(inn_clean, .keep_all = TRUE)

# Final coverage check
final_matched   <- our_inns[our_inns %in% ema_ma_dates_full$inn_clean]
final_unmatched <- our_inns[!our_inns %in% ema_ma_dates_full$inn_clean]

cat("=== FINAL COVERAGE ===\n")
cat("Matched:", length(final_matched), "/", length(our_inns), "\n")
if (length(final_unmatched) > 0) {
  cat("Still unmatched (will be NA in analysis):\n")
  print(final_unmatched)
  cat("These are intentionally left as NA — see ASSUMPTION [A3] in script header.\n")
}
cat("\n")

# Save the lookup table
write.csv(ema_ma_dates_full,
          here("data", "raw", "ema_ma_dates.csv"),
          row.names = FALSE, na = "")

cat("Saved: data/raw/ema_ma_dates.csv\n")
cat("Columns: inn_clean (join key), ema_ma_date (first EMA MA date), orphan_flag\n")
