# ============================================================
# combine_results.R
# Run after all cluster jobs complete
# ============================================================
library(dplyr); library(tidyr); library(ggplot2); library(patchwork)

# ── 1. Load all result files ──────────────────────────────────────────────────
files <- list.files("results", pattern = "\\.rds$",
                    recursive = TRUE, full.names = TRUE)
cat(sprintf("Found %d result files\n", length(files)))

results_list <- lapply(files, readRDS)
results      <- bind_rows(lapply(results_list, function(x) {
  # Drop the daily component for the main AR table
  x$daily <- NULL
  as.data.frame(x)
}))

cat(sprintf("Loaded %d rows (%d unique sims)\n",
            nrow(results),
            length(unique(paste(results$daycare_id, results$sim_id)))))

# ── 2. Compute baselines ──────────────────────────────────────────────────────
safe_prot <- function(ar, baseline) {
  ifelse(is.na(ar) | is.na(baseline) | baseline == 0, NA, 1 - ar / baseline)
}

baseline <- results %>%
  filter(vax_coverage == 0) %>%
  summarise(
    base_unvax_inf  = mean(ar_unvax_infant,  na.rm = TRUE),
    base_uvtod_dc   = mean(ar_unvax_tod_dc,  na.rm = TRUE),
    base_tod_hh     = mean(ar_toddler_hh,    na.rm = TRUE),
    base_adult      = mean(ar_adult,          na.rm = TRUE),
    base_elderly    = mean(ar_elderly,        na.rm = TRUE)
  )

# ── 3. Per-sim protection ─────────────────────────────────────────────────────
prot_per_sim <- results %>%
  filter(vax_coverage > 0) %>%
  mutate(
    vax_pct = factor(paste0(vax_coverage * 100, "%"),
                     levels = paste0(c(20,40,60,80,100), "%")),
    total_prot_infant   = 1 - ar_vax_infant   / baseline$base_unvax_inf,
    indirect_toddler_dc = 1 - ar_unvax_tod_dc / baseline$base_uvtod_dc,
    indirect_tod_hh_vaxhh = 1 - ar_tod_hh_vaxhh / baseline$base_tod_hh,
    indirect_adult      = 1 - ar_adult        / baseline$base_adult,
    indirect_elderly    = 1 - ar_elderly       / baseline$base_elderly
  )

# ── 4. Summarise to mean + 95% CI ─────────────────────────────────────────────
summarise_prot <- function(data, value_col) {
  data %>%
    group_by(vax_pct) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE) * 100,
      lo   = quantile(.data[[value_col]], 0.025, na.rm = TRUE) * 100,
      hi   = quantile(.data[[value_col]], 0.975, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(metric = value_col)
}

dc_metrics <- c("total_prot_infant", "indirect_toddler_dc")
hh_metrics <- c("indirect_tod_hh_vaxhh", "indirect_adult", "indirect_elderly")

dc_long <- bind_rows(lapply(dc_metrics, summarise_prot, data = prot_per_sim)) %>%
  mutate(metric = factor(metric, levels = dc_metrics,
                         labels = c("Direct: vaccinated infant",
                                    "Indirect: unvaccinated toddler (daycare)")))

hh_long <- bind_rows(lapply(hh_metrics, summarise_prot, data = prot_per_sim)) %>%
  mutate(metric = factor(metric, levels = hh_metrics,
                         labels = c("Toddler: HH has vaccinated infant",
                                    "Adult (parent)", "Elderly")))

# ── 5. Epidemic curves (Figures A/B) ─────────────────────────────────────────
# Load daily data only for curve files
curve_files <- list.files("results", pattern = "\\.rds$",
                          recursive = TRUE, full.names = TRUE)

daily_list <- lapply(curve_files, function(f) {
  x <- readRDS(f)
  if (is.null(x$daily)) return(NULL)
  x$daily$daycare_id   <- x$daycare_id
  x$daily$sim_id       <- x$sim_id
  x$daily$vax_coverage <- x$vax_coverage
  x$daily
})
daily_all <- bind_rows(Filter(Negate(is.null), daily_list))

curve_summary <- daily_all %>%
  group_by(vax_coverage, day) %>%
  summarise(
    dc_I_mean = mean(dc_I),
    dc_I_lo   = quantile(dc_I, 0.025),
    dc_I_hi   = quantile(dc_I, 0.975),
    hh_I_mean = mean(hh_I),
    hh_I_lo   = quantile(hh_I, 0.025),
    hh_I_hi   = quantile(hh_I, 0.975),
    .groups   = "drop"
  )

# ── 6. Save everything for plotting ──────────────────────────────────────────
saveRDS(list(
  prot_per_sim  = prot_per_sim,
  dc_long       = dc_long,
  hh_long       = hh_long,
  curve_summary = curve_summary,
  baseline      = baseline
), "results/combined_for_plots.rds")

cat("Saved: results/combined_for_plots.rds\n")
cat(sprintf("Total sims combined: %d across %d daycares\n",
            max(results$sim_id), length(unique(results$daycare_id))))

