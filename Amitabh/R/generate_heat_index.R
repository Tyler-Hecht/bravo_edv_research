#!/usr/bin/env Rscript
#
# generate_heat_index.R
#
# Main entry point. Generates daily, population-weighted heat index by HHS
# region, for any set of years, from NASA POWER.
#
# Pulls hourly temp/RH per county (power_hourly_pull.R), computes hourly
# heat index and takes the daily max (heat_index_formula.R) -- matches
# CDC's documented method. Then county -> state, plain mean (Eq. 2), and
# state -> region, population-weighted (Eq. 3 of the source manuscript).
#
# Scope: 50 states + DC (3,108 counties). Alaska/Hawaii excluded to match
# CDC's NLDAS-2-based coverage (contiguous US only) -- AK would need
# antimeridian-aware longitude handling, HI's county boundary includes
# uninhabited islands that would skew an area-weighted average. Neither is
# handled here.
#
# Cost: ~0.7s per county API call, ~35-40 min/year. Cached per county-year
# to INTERIM_CACHE_DIR, so a rerun on an already-pulled year is free.
#
# Requires: httr, jsonlite, dplyr, tidyr, lubridate, readr, sf

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(tidyr)
  library(lubridate); library(readr); library(sf)
})

# ============================================================================
# CONFIG -- edit these two lines to generate different years
# ============================================================================

YEARS <- c(2026)   # any year(s) NASA POWER has data for, e.g. c(2023, 2024, 2025)

OUTPUT_FILE <- "generated_heat_index.csv"   # written under OUTPUT_DIR below

# ---- paths (relative to this package's own layout; edit if you moved files) -

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  NULL
}
script_dir <- get_script_dir()
PACKAGE_ROOT <- if (!is.null(script_dir)) normalizePath(file.path(script_dir, "..")) else getwd()

source(file.path(PACKAGE_ROOT, "R", "heat_index_formula.R"))
source(file.path(PACKAGE_ROOT, "R", "power_hourly_pull.R"))

COUNTY_SHP        <- file.path(PACKAGE_ROOT, "inputs", "us_counties_shapefile", "cb_2018_us_county_500k.shp")
POP_FILE          <- file.path(PACKAGE_ROOT, "inputs", "state_population_2020_census.csv")
INTERIM_CACHE_DIR <- file.path(PACKAGE_ROOT, "cache", "power_hourly")   # per-county-year .rds cache
OUTPUT_DIR        <- file.path(PACKAGE_ROOT, "output")

# ============================================================================
# The 10 HHS regions and their states. AK (FIPS 02) and HI (FIPS 15) are
# deliberately absent -- see header note above.
# ============================================================================

get_region_states <- function() {
  tribble(
    ~region, ~state,                  ~state_abbr, ~statefp,
    1,       "Connecticut",           "CT",        "09",
    1,       "Maine",                 "ME",        "23",
    1,       "Massachusetts",         "MA",        "25",
    1,       "New Hampshire",         "NH",        "33",
    1,       "Rhode Island",          "RI",        "44",
    1,       "Vermont",               "VT",        "50",
    2,       "New Jersey",            "NJ",        "34",
    2,       "New York",              "NY",        "36",
    3,       "Delaware",              "DE",        "10",
    3,       "District of Columbia",  "DC",        "11",
    3,       "Maryland",              "MD",        "24",
    3,       "Pennsylvania",          "PA",        "42",
    3,       "Virginia",              "VA",        "51",
    3,       "West Virginia",         "WV",        "54",
    4,       "Alabama",               "AL",        "01",
    4,       "Florida",               "FL",        "12",
    4,       "Georgia",               "GA",        "13",
    4,       "Kentucky",              "KY",        "21",
    4,       "Mississippi",           "MS",        "28",
    4,       "North Carolina",        "NC",        "37",
    4,       "South Carolina",        "SC",        "45",
    4,       "Tennessee",             "TN",        "47",
    5,       "Illinois",              "IL",        "17",
    5,       "Indiana",               "IN",        "18",
    5,       "Michigan",              "MI",        "26",
    5,       "Minnesota",             "MN",        "27",
    5,       "Ohio",                  "OH",        "39",
    5,       "Wisconsin",             "WI",        "55",
    6,       "Arkansas",              "AR",        "05",
    6,       "Louisiana",             "LA",        "22",
    6,       "New Mexico",            "NM",        "35",
    6,       "Oklahoma",              "OK",        "40",
    6,       "Texas",                 "TX",        "48",
    7,       "Iowa",                  "IA",        "19",
    7,       "Kansas",                "KS",        "20",
    7,       "Missouri",              "MO",        "29",
    7,       "Nebraska",              "NE",        "31",
    8,       "Colorado",              "CO",        "08",
    8,       "Montana",               "MT",        "30",
    8,       "North Dakota",          "ND",        "38",
    8,       "South Dakota",          "SD",        "46",
    8,       "Utah",                  "UT",        "49",
    8,       "Wyoming",               "WY",        "56",
    9,       "Arizona",               "AZ",        "04",
    9,       "California",            "CA",        "06",
    9,       "Nevada",                "NV",        "32",
    10,      "Idaho",                 "ID",        "16",
    10,      "Oregon",                "OR",        "41",
    10,      "Washington",            "WA",        "53"
  )
}

# st_point_on_surface(), not a centroid -- a centroid can land outside a
# concave/coastal county (e.g. in a bay); point-on-surface always lands
# inside the polygon.
get_county_points <- function(shp_path, region_states) {
  counties_sf <- st_read(shp_path, quiet = TRUE) %>%
    filter(STATEFP %in% region_states$statefp) %>%
    st_transform(4326) %>%
    left_join(region_states, by = c("STATEFP" = "statefp"))

  pts <- st_point_on_surface(counties_sf)
  coords <- st_coordinates(pts)
  counties_sf %>%
    st_drop_geometry() %>%
    transmute(GEOID, region, state, state_abbr, lon = coords[, "X"], lat = coords[, "Y"])
}

# Generate one year of region-level heat index: pull/cache each county's
# hourly data, compute daily heat index, aggregate county -> state (Eq. 2)
# -> region (Eq. 3).
generate_year <- function(year, county_points, pop, cache_dir) {
  cat(sprintf("\n==== Year %d: %d counties ====\n", year, nrow(county_points)))

  t_start <- Sys.time()
  county_daily <- vector("list", nrow(county_points))
  failures <- character(0)

  for (i in seq_len(nrow(county_points))) {
    geoid <- county_points$GEOID[i]
    result <- tryCatch({
      hourly <- pull_county_hourly_cached(geoid, county_points$lat[i], county_points$lon[i], year, cache_dir)
      daily <- compute_daily_heat_index(hourly)
      daily %>% mutate(GEOID = geoid, region = county_points$region[i],
                        state = county_points$state[i], state_abbr = county_points$state_abbr[i])
    }, error = function(e) {
      message("  FAILED: ", geoid, " / ", year, " -- ", conditionMessage(e))
      failures <<- c(failures, geoid)
      NULL
    })
    county_daily[[i]] <- result
    if (i %% 200 == 0) {
      message(sprintf("  ... %d / %d counties done (%.1f min elapsed)",
                       i, nrow(county_points), as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
    }
  }

  cat(sprintf("Year %d: %.1f minutes, %d failures out of %d counties.\n",
              year, as.numeric(difftime(Sys.time(), t_start, units = "mins")), length(failures), nrow(county_points)))
  if (length(failures) > 0) {
    cat("Failed GEOIDs:", paste(failures, collapse = ", "), "\n")
    cat("These counties are simply missing from this year's output; rerun the script to retry them\n")
    cat("(already-successful counties are cached and will not be re-pulled).\n")
  }

  county_hi <- bind_rows(county_daily)

  # Eq. 2: county -> state, plain (unweighted) mean.
  state_hi <- county_hi %>%
    group_by(region, state, state_abbr, date) %>%
    summarise(hi_state_f = mean(heat_index_f, na.rm = TRUE), n_counties = n(), .groups = "drop")

  # Eq. 3: state -> region, population-weighted.
  state_hi %>%
    left_join(pop %>% select(state, population), by = "state") %>%
    group_by(region, date) %>%
    summarise(heat_index_f = sum(hi_state_f * population) / sum(population), .groups = "drop")
}

# ============================================================================
# run
# ============================================================================

region_states <- get_region_states()

pop <- read_csv(POP_FILE, show_col_types = FALSE, comment = "#")
missing_pop <- setdiff(region_states$state, pop$state)
if (length(missing_pop) > 0) stop("Missing population rows for: ", paste(missing_pop, collapse = ", "))

county_points <- get_county_points(COUNTY_SHP, region_states)

all_years <- bind_rows(lapply(YEARS, generate_year, county_points = county_points, pop = pop, cache_dir = INTERIM_CACHE_DIR))

output <- all_years %>%
  arrange(region, date) %>%
  transmute(`HHS Region` = region, Date = format(date, "%Y-%m-%d"), Region_Weighted_Heat_Index = heat_index_f)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(OUTPUT_DIR, OUTPUT_FILE)
write_csv(output, out_path)

cat(sprintf("\nDone. %d rows (%d regions x %d days) written to %s\n",
            nrow(output), length(unique(region_states$region)),
            nrow(output) / length(unique(region_states$region)), out_path))
