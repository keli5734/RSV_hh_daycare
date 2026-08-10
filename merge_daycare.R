# ============================================================================
# Merge per-centre packs into the combined per-coverage files that the
# plotting / summary scripts expect:  results1/RSV_rawsims_cov###.rds
# Run ONCE after all array tasks finish.
# ============================================================================
N_DAYCARES    <- 50L
VAX_LEVELS_MD <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
raw_dir <- file.path(getwd(), "results1", "raw")
out_dir <- file.path(getwd(), "results1")

for (vc in VAX_LEVELS_MD) {
  tag <- round(vc * 100)
  parts <- file.path(raw_dir, sprintf("RSV_rawsims_cov%03d_dc%03d.rds", tag, seq_len(N_DAYCARES)))
  miss  <- parts[!file.exists(parts)]
  if (length(miss) > 0)
    stop(sprintf("cov%03d: %d centre file(s) missing, e.g. %s -- did all array tasks finish?",
                 tag, length(miss), basename(miss[1])))
  combined <- do.call(c, lapply(parts, readRDS))     # concatenate lists of packs
  of <- file.path(out_dir, sprintf("RSV_rawsims_cov%03d.rds", tag))
  saveRDS(combined, of)
  cat(sprintf("cov %3.0f%%: merged %d centres -> %d packs -> %s\n",
              vc * 100, N_DAYCARES, length(combined), basename(of)))
}
cat("\nMerge complete. results1/RSV_rawsims_cov*.rds ready for plotting/summary scripts.\n")
cat("(You may delete results1/raw/ to reclaim space once verified.)\n")
