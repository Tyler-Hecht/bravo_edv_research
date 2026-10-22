#!/usr/bin/env Rscript
#
# heat_index_formula.R
#
# NWS "modified Rothfusz" heat index. Core regression is Rothfusz (1990)
# (Eq. 1 in the source manuscript); "modified" = NWS's standard branches on
# top of it: a linear estimate for mild conditions (the raw regression only
# holds up in hot/humid ranges), plus small corrections at very low/high
# humidity. Matches the published weather.gov algorithm -- not diffed
# line-by-line against CDC's own source, but tracks CDC's published values
# closely in testing. If a future comparison drifts, check the branch
# thresholds/corrections below first.
#
# Only place in the pipeline that does this math -- everything else calls
# compute_heat_index().

# t_f: air temp (F). rh_pct: relative humidity, 0-100. Vectorized.
compute_heat_index <- function(t_f, rh_pct) {
  simple <- 0.5 * (t_f + 61.0 + (t_f - 68.0) * 1.2 + rh_pct * 0.094)

  # NWS switches on the average of `simple` and t_f, not `simple` alone --
  # using simple alone misclassifies some mild-humidity days into the full
  # regression.
  avg_est <- (simple + t_f) / 2

  full <- -42.379 +
    2.04901523 * t_f +
    10.14333127 * rh_pct -
    0.22475541 * t_f * rh_pct -
    6.83783e-3 * t_f^2 -
    5.481717e-2 * rh_pct^2 +
    1.22874e-3 * t_f^2 * rh_pct +
    8.5282e-4 * t_f * rh_pct^2 -
    1.99e-6 * t_f^2 * rh_pct^2

  # Low-RH (hot+dry) and high-RH (hot+very humid) corrections. The two
  # trigger windows don't overlap, so sequential ifelse() calls are fine.
  # pmax(0, ...) just guards sqrt() outside low_rh_adj's valid window --
  # use_full gates it to 80<=t_f<=112 below regardless.
  low_rh_adj  <- ((13 - rh_pct) / 4) * sqrt(pmax(0, (17 - abs(t_f - 95)) / 17))
  high_rh_adj <- ((rh_pct - 85) / 10) * ((87 - t_f) / 5)

  use_full <- avg_est >= 80
  full_adj <- full
  full_adj <- ifelse(use_full & rh_pct < 13 & t_f >= 80 & t_f <= 112, full - low_rh_adj, full_adj)
  full_adj <- ifelse(use_full & rh_pct > 85 & t_f >= 80 & t_f <= 87, full + high_rh_adj, full_adj)

  ifelse(t_f < 40, t_f, ifelse(use_full, full_adj, simple))
}
