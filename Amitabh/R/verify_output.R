#!/usr/bin/env Rscript
#
# verify_output.R
#
# Sanity-checks a file produced by generate_heat_index.R:
#   1. Completeness (every region has every expected day) and bounds (no
#      negative/implausible values) -- always run.
#   2. Accuracy vs. CDC's real data -- only for years CDC has published
#      (through 2022). Set CDC_REFERENCE_FILE to run it; leave it NA for
#      years with no ground truth (the normal case, since that's the whole
#      reason this pipeline exists). Confidence in an unpublished year
#      comes from this check passing on a year CDC HAS published, not from
#      anything checkable directly.
#
# Edit the two paths below, then run.

suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(readr)
})

# ============================================================================
# CONFIG -- edit these
# ============================================================================

GENERATED_FILE <- "output/generated_heat_index.csv"   # a file produced by generate_heat_index.R
CDC_REFERENCE_FILE <- NA   # set to a CDC file's path if the generated year(s) overlap it; else NA

BOUNDS_MIN <- 0     # heat index below this is flagged as implausible
BOUNDS_MAX <- 130   # heat index above this is flagged as implausible

# ============================================================================
# Check 1: completeness -- every region, every expected day
# ============================================================================

# Every region should have 365 days (366 in a leap year) for each year present.
check_completeness <- function(d) {
  d %>%
    mutate(yr = year(date)) %>%
    group_by(region, yr) %>%
    summarise(actual_days = n(), .groups = "drop") %>%
    mutate(
      expected_days = if_else(leap_year(yr), 366L, 365L),
      flag = if_else(actual_days != expected_days, "FLAG: day count mismatch", "ok")
    )
}

check_bounds <- function(d, lo = BOUNDS_MIN, hi = BOUNDS_MAX) {
  d %>%
    group_by(region) %>%
    summarise(
      n_days = n(), min_hi = round(min(heat_index_f), 2), max_hi = round(max(heat_index_f), 2),
      n_below_lo = sum(heat_index_f < lo), n_above_hi = sum(heat_index_f > hi),
      .groups = "drop"
    ) %>%
    mutate(flag = if_else(n_below_lo > 0 | n_above_hi > 0, "FLAG: implausible value(s)", "ok"))
}

# Bias/RMSE/MAE/correlation vs. CDC data, per region, over overlapping dates.
check_against_cdc <- function(generated, cdc_file) {
  cdc <- read_csv(cdc_file, show_col_types = FALSE) %>%
    rename(region = `HHS Region`, hi_cdc = Region_Weighted_Heat_Index) %>%
    mutate(date = as_date(Date))

  comp <- generated %>%
    inner_join(cdc %>% select(region, date, hi_cdc), by = c("region", "date")) %>%
    mutate(diff = heat_index_f - hi_cdc)

  if (nrow(comp) == 0) {
    message("No overlapping region-dates between the generated file and the CDC reference file -- nothing to compare.")
    return(invisible(NULL))
  }

  comp %>%
    group_by(region) %>%
    summarise(
      days_compared = n(),
      bias = round(mean(diff), 2), rmse = round(sqrt(mean(diff^2)), 2), mae = round(mean(abs(diff)), 2),
      correlation = round(cor(heat_index_f, hi_cdc), 4),
      .groups = "drop"
    )
}

# ============================================================================
# run
# ============================================================================

d <- read_csv(GENERATED_FILE, show_col_types = FALSE) %>%
  rename(region = `HHS Region`, heat_index_f = Region_Weighted_Heat_Index) %>%
  mutate(date = as_date(Date))

cat("==== Completeness ====\n")
completeness <- check_completeness(d)
print(completeness, n = 30)
cat(sprintf("\n%d flag(s)\n", sum(completeness$flag != "ok")))

cat("\n==== Bounds ====\n")
bounds <- check_bounds(d)
print(bounds, n = 15)
cat(sprintf("\n%d flag(s)\n", sum(bounds$flag != "ok")))

if (!is.na(CDC_REFERENCE_FILE)) {
  cat("\n==== Accuracy vs. CDC (only meaningful where both cover the same region-dates) ====\n")
  accuracy <- check_against_cdc(d, CDC_REFERENCE_FILE)
  if (!is.null(accuracy)) print(accuracy, n = 15)
} else {
  cat("\nCDC_REFERENCE_FILE not set -- skipping the accuracy check (expected if the generated year(s)\n")
  cat("have no CDC ground truth, which is the normal case for this pipeline).\n")
}
