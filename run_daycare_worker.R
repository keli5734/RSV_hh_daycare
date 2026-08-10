# ============================================================================
# Array worker: run ONE daycare centre (all coverages, all sims) and save its
# per-centre raw packs. Centre index comes from $SLURM_ARRAY_TASK_ID (1-based).
# Reuses rsv_multidaycare_rawdata.R's seeding so output == serial run.
#   Output: results1/raw/RSV_rawsims_cov###_dc###.rds   (one per coverage)
# ============================================================================
suppressWarnings(suppressMessages({
  RUN_ANALYSIS <- FALSE
  source("rsv_model_M.R")
}))

## ---- config (match your serial run) ---------------------------------------
N_DAYCARES <- 50L
params$n_daycares <- N_DAYCARES
params$n_sims     <- 100L

VAX_LEVELS_MD <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
POP_SEED_BASE <- params$pop_seed     # 20261234
DYN_SEED_BASE <- 7000000L            # same as driver
RAW_COMPONENTS <- c("daily","agents_final","infection_log")

## ---- which centre am I? ----------------------------------------------------
d <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
if (is.na(d) || d < 1L || d > N_DAYCARES)
  stop(sprintf("SLURM_ARRAY_TASK_ID=%s out of range 1..%d", Sys.getenv("SLURM_ARRAY_TASK_ID"), N_DAYCARES))

raw_dir <- file.path(getwd(), "results1", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

pack_raw <- function(res, d, s, vc, components) {
  out <- list(daycare_id = d, sim = s, vax = vc)
  if ("daily"         %in% components) out$daily         <- res$daily
  if ("agents_final"  %in% components) out$agents_final  <- res$agents_final
  if ("infection_log" %in% components) out$infection_log <- res$infection_log
  out
}

cat(sprintf("=== centre %d / %d : %d sims x %d coverages ===\n",
            d, N_DAYCARES, params$n_sims, length(VAX_LEVELS_MD)))
t0 <- Sys.time()

p <- params; p$pop_seed <- POP_SEED_BASE + d        # centre-specific population
for (vc in VAX_LEVELS_MD) {
  cov_off <- as.integer(round(vc * 1000)) * 100000L
  raw <- vector("list", params$n_sims)
  for (s in seq_len(params$n_sims)) {
    set.seed(DYN_SEED_BASE + cov_off + d * 10000L + s)   # IDENTICAL to serial driver
    res <- run_simulation(p, vax_coverage = vc)
    raw[[s]] <- pack_raw(res, d, s, vc, RAW_COMPONENTS)
  }
  rf <- file.path(raw_dir, sprintf("RSV_rawsims_cov%03d_dc%03d.rds", round(vc * 100), d))
  saveRDS(raw, rf)
  cat(sprintf("  cov %3.0f%% -> %s\n", vc * 100, basename(rf)))
  rm(raw); gc(verbose = FALSE)
}
cat(sprintf("=== centre %d done in %.1f min ===\n", d, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
