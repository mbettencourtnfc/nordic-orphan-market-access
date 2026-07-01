# =============================================================================
# 00_fetch_ema_ma_dates.R
# Fetch EMA first-marketing-authorisation dates for our 51 drugs
# =============================================================================
# Downloads the EMA authorised medicines Excel (updated nightly by EMA),
# extracts INN + first authorisation date, and saves a lookup CSV.
#
# Run this ONCE before running 01_data_cleaning.R.
# Output: data/raw/ema_ma_dates.csv
# =============================================================================

library(readxl)
library(dplyr)
library(stringr)
library(here)
library(lubridate)

# -----------------------------------------------------------------------------
# 1. Download the EMA medicines Excel
# -----------------------------------------------------------------------------
ema_url  <- "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx"
dest_file <- here("data", "raw", "ema_medicines_raw.xlsx")

cat("Downloading EMA medicines data...\n")
download.file(ema_url, destfile = dest_file, mode = "wb", quiet = FALSE)
cat("Download complete.\n\n")

# -----------------------------------------------------------------------------
# 2. Read and extract INN + first marketing authorisation date
# -----------------------------------------------------------------------------
ema_raw <- read_excel(dest_file, skip = 8)

ema_ma_dates <- ema_raw |>
  filter(Category == "Human", `Medicine status` == "Authorised") |>
  select(
    inn_ema          = `International non-proprietary name (INN) / common name`,
    ma_date_raw      = `Marketing authorisation date`,
    orphan_flag      = `Orphan medicine`
  ) |>
  filter(!is.na(inn_ema), !is.na(ma_date_raw)) |>
  mutate(
    inn_clean       = str_to_lower(str_trim(inn_ema)),
    # Column is numeric (Excel serial date) or character depending on export
    ema_ma_date     = case_when(
      is.numeric(ma_date_raw) ~ as_date(as.numeric(ma_date_raw), origin = "1899-12-30"),
      TRUE                    ~ dmy(as.character(ma_date_raw))
    )
  ) |>
  filter(!is.na(ema_ma_date)) |>
  # A drug can have multiple authorisations (extensions); keep the earliest
  group_by(inn_clean) |>
  slice_min(ema_ma_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(inn_clean, ema_ma_date, orphan_flag)

cat("EMA MA dates extracted:", nrow(ema_ma_dates), "unique INNs\n")
cat("Date range:", format(min(ema_ma_dates$ema_ma_date)), "to",
    format(max(ema_ma_dates$ema_ma_date)), "\n\n")

# Check how many of our 51 drugs match
master <- read_excel(here("data", "raw", "nordic_orphan_master.xlsx"), sheet = "Master Data")
our_inns <- str_to_lower(str_trim(unique(master$INN)))
matched  <- our_inns[our_inns %in% ema_ma_dates$inn_clean]
cat("Drugs matched from our dataset:", length(matched), "/", length(our_inns), "\n")
if (length(our_inns) > length(matched)) {
  cat("Unmatched INNs:\n")
  print(our_inns[!our_inns %in% ema_ma_dates$inn_clean])
}
cat("\n")

# Sample of matched dates
cat("Sample matched dates:\n")
master |>
  mutate(inn_clean = str_to_lower(str_trim(INN))) |>
  distinct(Drug, inn_clean) |>
  left_join(ema_ma_dates, by = "inn_clean") |>
  filter(!is.na(ema_ma_date)) |>
  head(10) |>
  print()
cat("\n")

# -----------------------------------------------------------------------------
# 3. Fuzzy-match unmatched drugs
# -----------------------------------------------------------------------------
unmatched_inns <- our_inns[!our_inns %in% ema_ma_dates$inn_clean]

# (a) Combination drugs: match on first component before "/"
combo_fix <- data.frame(inn_master = unmatched_inns[grepl("/", unmatched_inns)]) |>
  mutate(
    first_component = str_to_lower(str_trim(str_extract(inn_master, "^[^/]+")))
  ) |>
  left_join(ema_ma_dates |> select(first_component = inn_clean, ema_ma_date, orphan_flag),
            by = "first_component") |>
  filter(!is.na(ema_ma_date)) |>
  select(inn_clean = inn_master, ema_ma_date, orphan_flag)

cat("Combination drug matches found:", nrow(combo_fix), "\n")
print(combo_fix)
cat("\n")

# (b) Partial string match for gene therapies (e.g. "betibeglogene" in EMA INN)
still_unmatched <- unmatched_inns[!unmatched_inns %in% combo_fix$inn_clean]
still_unmatched <- still_unmatched[!grepl("/", still_unmatched)]

fuzzy_fix <- lapply(still_unmatched, function(our_inn) {
  # Take first word of our INN (most distinctive part)
  first_word <- str_extract(our_inn, "^\\S+")
  hit <- ema_ma_dates |>
    filter(str_detect(inn_clean, fixed(first_word, ignore_case = TRUE))) |>
    slice_min(ema_ma_date, n = 1, with_ties = FALSE)
  if (nrow(hit) == 0) return(NULL)
  tibble(inn_clean = our_inn, ema_ma_date = hit$ema_ma_date, orphan_flag = hit$orphan_flag,
         matched_ema_inn = hit$inn_clean)
}) |> bind_rows()

cat("Fuzzy gene therapy matches found:", nrow(fuzzy_fix), "\n")
if (nrow(fuzzy_fix) > 0) print(fuzzy_fix)
cat("\n")

# (c) Manual overrides — verified from individual EMA medicine pages
#     Includes drugs that: (i) use slash-combination INN format, (ii) had
#     conditional/withdrawn MA (so filtered out above), or (iii) are
#     listed under a slightly different INN in the EMA file.
#     Source URL pattern: ema.europa.eu/en/medicines/human/EPAR/<medicine-name>
manual_overrides <- tribble(
  ~inn_clean,                              ~ema_ma_date,           ~orphan_flag, ~note,
  "lumacaftor/ivacaftor",                  as_date("2015-11-19"),  "Yes",        "Orkambi — EMA conditional MA",
  "elexacaftor/tezacaftor/ivacaftor",      as_date("2020-08-21"),  "Yes",        "Kaftrio",
  "betibeglogene autotemcel",              as_date("2019-08-29"),  "Yes",        "Zynteglo — conditional MA (later withdrawn)",
  "atidarsagene autotemcel",               as_date("2020-12-16"),  "Yes",        "Libmeldy",
  "elivaldogene autotemcel",               as_date("2021-07-26"),  "Yes",        "Skysona",
  "fidanacogene elaparvovec",              as_date("2024-09-19"),  "Yes",        "Beqvez",
  "ataluren",                              as_date("2014-07-31"),  "Yes",        "Translarna — conditional MA (renewal refused 2023)",
  # delandistrogene moxeparvovec (Elevidys): FDA 2023; EMA review ongoing as of 2024 — leave NA
  # birch triterpenes / bjørkeneverekstrakt: Norwegian-only product, no EMA central authorisation — leave NA
) |> select(inn_clean, ema_ma_date, orphan_flag)

cat("Manual overrides added:", nrow(manual_overrides), "\n\n")

# -----------------------------------------------------------------------------
# 4. Merge all matches and save
# -----------------------------------------------------------------------------
ema_ma_dates_full <- bind_rows(
  ema_ma_dates |> select(inn_clean, ema_ma_date, orphan_flag),
  combo_fix,
  fuzzy_fix |> select(inn_clean, ema_ma_date, orphan_flag),
  manual_overrides
) |>
  distinct(inn_clean, .keep_all = TRUE)

# Final coverage check
final_matched <- our_inns[our_inns %in% ema_ma_dates_full$inn_clean]
cat("Final coverage:", length(final_matched), "/", length(our_inns), "drugs\n\n")

write.csv(ema_ma_dates_full, here("data", "raw", "ema_ma_dates.csv"),
          row.names = FALSE, na = "")
cat("Saved: data/raw/ema_ma_dates.csv\n")
