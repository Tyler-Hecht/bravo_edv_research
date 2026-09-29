#!/usr/bin/env Rscript
#
# power_hourly_pull.R
#
# Pulls hourly temp/RH from NASA POWER and collapses one county-year down
# to daily heat index, same procedure CDC uses: heat index computed per
# hour, day's value = max across that day's hours.
#
# Point-based, not area-averaged: NASA POWER's hourly endpoint only takes
# single points, no bounding box (unlike its daily endpoint). Each county
# is one representative point -- see st_point_on_surface() in
# generate_heat_index.R. Disclosed simplification, see README's "Known
# limitations".
#
# Every county-year pull is cached to its own .rds, so an interrupted run
# resumes for free.

# One NASA POWER hourly API call. start/end as "YYYYMMDD". Returns raw
# t2m_c/rh_pct, fill values (-999) left as-is -- cleaned in the caller.
fetch_power_hourly_point <- function(lat, lon, start, end) {
  url <- "https://power.larc.nasa.gov/api/temporal/hourly/point"
  resp <- httr::RETRY(
    "GET", url,
    query = list(
      parameters = "T2M,RH2M", community = "AG",
      longitude = lon, latitude = lat, start = start, end = end, format = "JSON"
    ),
    times = 5, pause_base = 2, pause_cap = 20,
    terminate_on = c(400, 401, 403, 404)
  )
  httr::stop_for_status(resp)
  parsed <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  p <- parsed$properties$parameter
  tibble::tibble(
    datetime = names(p$T2M),
    t2m_c = unlist(p$T2M),
    rh_pct = unlist(p$RH2M[names(p$T2M)])
  )
}

# Pull (or load from cache) one county's full hourly year -- one API call
# covers Jan 1 - Dec 31, so one call per county-year.
pull_county_hourly_cached <- function(geoid, lat, lon, year, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, sprintf("%s_%d.rds", geoid, year))

  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  message("  pulling hourly data: county ", geoid, ", ", year, " ...")
  start <- sprintf("%d0101", year)
  end   <- sprintf("%d1231", year)
  raw <- fetch_power_hourly_point(lat, lon, start, end)

  clean <- dplyr::mutate(
    raw,
    t2m_c = ifelse(t2m_c <= -990, NA_real_, t2m_c),   # NASA POWER fill value
    rh_pct = ifelse(rh_pct <= -990, NA_real_, rh_pct),
    date = as.Date(substr(datetime, 1, 8), format = "%Y%m%d"),
    t2m_f = t2m_c * 9 / 5 + 32
  )
  clean <- dplyr::filter(clean, !is.na(t2m_f), !is.na(rh_pct))

  saveRDS(clean, cache_file)
  clean
}

# Collapse cleaned hourly data to daily heat index (max across the day's
# hours). A day missing some hours just takes the max over whatever's
# present rather than warning -- verify_output.R's completeness check
# catches that case.
compute_daily_heat_index <- function(hourly_df) {
  hourly_df %>%
    dplyr::mutate(hi_hourly_f = compute_heat_index(t2m_f, rh_pct)) %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(heat_index_f = max(hi_hourly_f), .groups = "drop")
}
