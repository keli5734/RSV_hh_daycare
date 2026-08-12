# ============================================================================
# reconstruct_contact_matrix_M.R
#
# Aggregates 17x17 Epistorm-Mix MEAN-number-of-contacts (M) matrices
# (5-year bands: 0-4, 5-9, ..., 75-79, 80+) into three sub-matrices used by
# the RSV daycare-household model:
#
#   M_infant_room  : infant (<1)  x toddler (1-4)  x staff (18-64)
#   M_toddler_room : toddler (1-4) x staff (18-64)
#   M_household    : infant (<1)  x toddler (1-4)  x adult (18-64) x elderly (65+)
#
# M[i,j] = mean number of contacts that a participant in age band i reports
#          with individuals in age band j (per day, setting-specific).
#
# WHY M (not F): the transmission model builds its force of infection as
#   lambda_i = sum_j  C[i,j] * (I_j / N_j)
# i.e. (mean contacts with j) x (prevalence in j). That structure expects a
# COUNT in C[i,j]; the M-matrix supplies exactly that. (F-matrices are a
# per-capita contact PROBABILITY and do not belong in that I_j/N_j form.)
#
# KEY AGGREGATION RULE for M-matrices, merging bands -> super-groups I, J:
#   - contact columns  -> SUM:   contacts with the merged group J = sum_{j in J} M[i,j]
#   - participant rows -> population-weighted MEAN over the bands in I
#   =>  M[I,J] = sum_{i in I} (N_i / N_I) * sum_{j in J} M[i,j]
#   In R:  sum(pop_rows * rowSums(sub)) / sum(pop_rows)
#   (Population enters ONLY on the rows. Columns are a plain sum.)
# ============================================================================

# ── 0. USER INPUT ─────────────────────────────────────────────────────────────
# NOTE: confirm the exact M-matrix filenames in your Epistorm-Mix download.
#       They mirror the F naming, with "M" in place of "F".
FILE_PATH  <- "contact_matrix/School-M-by5_80-matrix.csv"     # daycare proxy
FILE_PATH2 <- "contact_matrix/Household-M-by5_80-matrix.csv"  # household

VERBOSE             <- TRUE   # print matrices + diagnostics
BALANCE_RECIPROCITY <- TRUE  # if TRUE, enforce N_i*M[i,j] = N_j*M[j,i] (see §10b)

# ── 1. BAND DEFINITIONS ───────────────────────────────────────────────────────
# Row i = respondent (participant) age band; Column j = contact age band.
BAND_LABELS <- c(
  "0-4",  "5-9",  "10-14","15-19",
  "20-24","25-29","30-34","35-39",
  "40-44","45-49","50-54","55-59",
  "60-64","65-69","70-74","75-79","80+"
)
N_BANDS <- length(BAND_LABELS)   # 17

# ── 2. US POPULATION BY 5-YEAR BAND ───────────────────────────────────────────
# US Census Bureau, National Population Estimates 2022 (thousands).
# Only the *ratios* matter for the row-weighting. Update vintage as preferred.
POP_BY_BAND <- c(
  "0-4"   = 19500, "5-9"   = 20200, "10-14" = 21000, "15-19" = 21500,
  "20-24" = 21800, "25-29" = 22800, "30-34" = 23000, "35-39" = 22200,
  "40-44" = 21000, "45-49" = 20200, "50-54" = 21800, "55-59" = 22300,
  "60-64" = 21200, "65-69" = 18500, "70-74" = 15500, "75-79" = 9900,
  "80+"   = 13100
)
stopifnot(identical(names(POP_BY_BAND), BAND_LABELS))   # guard against typos

# ── 3. MAP BROAD ROLES TO BAND INDICES ────────────────────────────────────────
#  infant  : ages  0 – <1  -> inside band 1 ("0-4"), split internally (col only)
#  toddler : ages  1 –  4  -> inside band 1 ("0-4"), split internally (col only)
#  staff   : ages 18 – 64  -> bands 5-13 ("20-24" … "60-64")
#  elderly : ages 65+      -> bands 14-17 ("65-69","70-74","75-79","80+")
#  Bands 2-4 ("5-9","10-14","15-19") = school-age children, not in the model.
ROLE_BANDS <- list(
  infant  = 1L,
  toddler = 1L,
  staff   = 5:13,
  elderly = 14:17
)

# ── 4. WITHIN-BAND SPLIT WEIGHTS: infant vs toddler inside "0-4" ──────────────
# "0-4" spans 5 single-year ages (0,1,2,3,4). Under a uniform within-band age
# distribution: infant (age 0) = 1/5 of the band, toddler (ages 1-4) = 4/5.
#
# For M (a COUNT), the split is a clean apportionment:
#   * As a CONTACT (column): if a band-1 person makes M[i,"0-4"] contacts with
#     band 1, and contacts are proportional to within-band sub-population, then
#       contacts with infant  = M[i,"0-4"] * 1/5
#       contacts with toddler = M[i,"0-4"] * 4/5
#   * As a PARTICIPANT (row): the data only give the band average. With no
#     finer resolution we assume an infant participant and a toddler participant
#     each make the band-average number of contacts (row rate unchanged).
#
# CONSEQUENCE (documented limitation): because the row rate is shared, the
# infant and toddler ROWS of each output matrix come out identical. The model
# differentiates infants from toddlers downstream via phi (susceptibility) and
# kappa (infectivity), not via the contact matrix. Replace 1/5, 4/5 with
# single-year census shares if you want a non-uniform split.
w_infant  <- 1 / 5
w_toddler <- 4 / 5

# ── 5. LOAD A RAW 17x17 M-MATRIX ──────────────────────────────────────────────
# Excel autoformatting turns some band labels into dates when the CSV is opened
# and saved in Excel: "5-9" -> "9-May" (or "5-Sep"), "10-14" -> "14-Oct"
# (or "10-Oct"). The numeric VALUES and column ORDER are unaffected — only the
# text labels. This restores them so the order-safety check still works.
# (The loader reassigns BAND_LABELS by position regardless, so results are
#  correct either way; this just keeps the safety net meaningful.)
normalize_band_labels <- function(x) {
  x <- trimws(as.character(x))
  recode <- c(
    "9-May"  = "5-9",   "5-Sep"  = "5-9",   "May-09" = "5-9",   "Sep-05" = "5-9",
    "14-Oct" = "10-14", "10-Oct" = "10-14", "Oct-14" = "10-14", "Oct-10" = "10-14"
  )
  ifelse(x %in% names(recode), unname(recode[x]), x)
}

load_m_matrix <- function(path) {
  raw <- read.csv(path, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)

  # Layout detection (N_BANDS = 17 for the by5_80 matrices):
  #   (a) Wide N x (N+1): first column is a row-label, remaining N are contact bands
  #   (b) Wide N x N: all numeric, no label column
  #   (c) Long: columns resp_age / cont_age / value
  file_labels <- NULL
  col_labels  <- NULL
  if (ncol(raw) == N_BANDS + 1 || ncol(raw) == N_BANDS) {
    if (is.character(raw[, 1])) {
      file_labels <- normalize_band_labels(raw[, 1])
      col_labels  <- normalize_band_labels(names(raw)[-1])  # header, minus label col
      mat <- as.matrix(raw[, -1])
    } else {
      col_labels <- normalize_band_labels(names(raw))
      mat <- as.matrix(raw)
    }
  } else if (all(c("resp_age", "cont_age", "value") %in% tolower(names(raw)))) {
    if (!requireNamespace("tidyr", quietly = TRUE))
      stop("Long-format CSV detected but the 'tidyr' package is not installed.")
    names(raw) <- tolower(names(raw))
    wide <- tidyr::pivot_wider(raw, names_from = cont_age, values_from = value)
    file_labels <- normalize_band_labels(wide[[1]])
    col_labels  <- normalize_band_labels(names(wide)[-1])
    mat <- as.matrix(wide[, -1])
  } else {
    stop(
      "Cannot recognise CSV layout. Expected an N x N (wide) matrix or a ",
      "long-format table with columns resp_age / cont_age / value.\n",
      "Columns found: ", paste(names(raw), collapse = ", ")
    )
  }

  mat <- apply(mat, 2, as.numeric)   # strips dimnames; reassigned below

  if (nrow(mat) != N_BANDS || ncol(mat) != N_BANDS)
    stop("Matrix must be ", N_BANDS, "x", N_BANDS,
         " after loading. Got: ", nrow(mat), "x", ncol(mat))

  # Order check (after un-mangling Excel date artifacts). Rows AND columns must
  # match the canonical 0-4 … 80+ order; warn loudly if not, since the loader
  # maps by position and a genuinely reordered file would be silently wrong.
  if (!is.null(file_labels) && length(file_labels) == N_BANDS &&
      !identical(file_labels, BAND_LABELS)) {
    warning(
      "File ROW labels differ from expected band order — verify ordering.\n",
      "  file:     ", paste(file_labels, collapse = ", "), "\n",
      "  expected: ", paste(BAND_LABELS, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(col_labels) && length(col_labels) == N_BANDS &&
      !identical(col_labels, BAND_LABELS)) {
    warning(
      "File COLUMN labels differ from expected band order — verify ordering.\n",
      "  file:     ", paste(col_labels, collapse = ", "), "\n",
      "  expected: ", paste(BAND_LABELS, collapse = ", "),
      call. = FALSE
    )
  }

  rownames(mat) <- BAND_LABELS
  colnames(mat) <- BAND_LABELS
  mat
}

M_raw  <- load_m_matrix(FILE_PATH)    # school / daycare proxy
M_raw2 <- load_m_matrix(FILE_PATH2)   # household

if (VERBOSE) {
  cat("\n── Raw 17x17 M-matrix, school/daycare (first 5 rows/cols) ──\n")
  print(round(M_raw[1:5, 1:5], 4))
}

# ── 6. CORE AGGREGATION (M semantics): sum cols, pop-weight rows ──────────────
aggregate_M <- function(M, role_row, role_col, pop = POP_BY_BAND) {
  band_rows <- ROLE_BANDS[[role_row]]
  band_cols <- ROLE_BANDS[[role_col]]
  pop_rows  <- pop[band_rows]

  sub <- M[band_rows, band_cols, drop = FALSE]

  # contacts with merged contact group  = SUM over contact columns
  # participant aggregation             = population-weighted MEAN over rows
  base <- sum(pop_rows * rowSums(sub)) / sum(pop_rows)

  # within-band CONTACT (column) split for infant / toddler (band 1 only)
  w_col <- if (role_col == "infant")  w_infant
  else     if (role_col == "toddler") w_toddler
  else      1

  base * w_col
}

build_sub_matrix <- function(M, roles, pop = POP_BY_BAND) {
  n <- length(roles)
  out <- matrix(0, n, n, dimnames = list(roles, roles))
  for (i in seq_len(n))
    for (j in seq_len(n))
      out[i, j] <- aggregate_M(M, roles[i], roles[j], pop)
  out
}

# ── 7. ROLE POPULATIONS (for reciprocity diagnostics / balancing) ─────────────
role_pop <- c(
  infant  = unname(w_infant  * POP_BY_BAND["0-4"]),
  toddler = unname(w_toddler * POP_BY_BAND["0-4"]),
  staff   = unname(sum(POP_BY_BAND[ROLE_BANDS$staff])),
  adult   = unname(sum(POP_BY_BAND[ROLE_BANDS$staff])),   # adult == staff age range
  elderly = unname(sum(POP_BY_BAND[ROLE_BANDS$elderly]))
)

# Reciprocity-balance an aggregated matrix so that N_i*M[i,j] = N_j*M[j,i].
# Empirical M-matrices are only approximately reciprocal; balancing removes
# residual asymmetry from smoothing/sampling. This is a modeling choice.
balance_reciprocity <- function(M, npop) {
  roles <- rownames(M)
  np <- npop[roles]
  Mb <- M
  for (a in roles) for (b in roles)
    Mb[a, b] <- (np[a] * M[a, b] + np[b] * M[b, a]) / (2 * np[a])
  Mb
}

# ── 8. CONSTRUCT THE THREE TARGET MATRICES ────────────────────────────────────
M_infant_room  <- build_sub_matrix(M_raw,  c("infant", "toddler", "staff"))
M_toddler_room <- build_sub_matrix(M_raw,  c("toddler", "staff"))
M_household    <- build_sub_matrix(M_raw2, c("infant", "toddler", "staff", "elderly"))
rownames(M_household) <- c("infant", "toddler", "adult", "elderly")
colnames(M_household) <- c("infant", "toddler", "adult", "elderly")

if (BALANCE_RECIPROCITY) {
  M_infant_room  <- balance_reciprocity(M_infant_room,  role_pop)
  M_toddler_room <- balance_reciprocity(M_toddler_room, role_pop)
  M_household    <- balance_reciprocity(M_household,    role_pop)
}

# ── 9. PRINT & INSPECT ────────────────────────────────────────────────────────
if (VERBOSE) {
  div <- strrep("=", 56)
  cat("\n", div, "\n M_infant_room  (infant x toddler x staff)\n", div, "\n", sep = "")
  print(round(M_infant_room, 5))
  cat("\n", div, "\n M_toddler_room  (toddler x staff)\n", div, "\n", sep = "")
  print(round(M_toddler_room, 5))
  cat("\n", div, "\n M_household  (infant x toddler x adult x elderly)\n", div, "\n", sep = "")
  print(round(M_household, 5))
}

# ── 10a. AGGREGATION-FIX DIAGNOSTIC: column SUM vs column MEAN ─────────────────
# The previous (F-style) script averaged contact columns, understating any
# multi-band contact group by ~(number of bands it spans): staff ~9x, elderly
# ~3x. This shows the magnitude of that fix on the current M data.
if (VERBOSE) {
  aggregate_M_avgcols <- function(M, role_row, role_col, pop = POP_BY_BAND) {
    br <- ROLE_BANDS[[role_row]]; bc <- ROLE_BANDS[[role_col]]
    pr <- pop[br]
    sub <- M[br, bc, drop = FALSE]
    base <- sum(pr * rowMeans(sub)) / sum(pr)   # MEAN over cols (old behavior)
    w_col <- if (role_col == "infant") w_infant else
      if (role_col == "toddler") w_toddler else 1
    base * w_col
  }
  build_avgcols <- function(M, roles, pop = POP_BY_BAND) {
    n <- length(roles); out <- matrix(0, n, n, dimnames = list(roles, roles))
    for (i in seq_len(n)) for (j in seq_len(n))
      out[i, j] <- aggregate_M_avgcols(M, roles[i], roles[j], pop)
    out
  }
  ir_avg <- build_avgcols(M_raw, c("infant", "toddler", "staff"))
  cat("\n── Ratio (correct SUM-cols / old MEAN-cols), M_infant_room ──\n")
  cat("   (≈ number of contact bands summed; staff col should be ~9x)\n")
  print(round(M_infant_room / ir_avg, 3))
}

# ── 10b. RECIPROCITY DIAGNOSTIC ───────────────────────────────────────────────
reciprocity_report <- function(M, npop, name) {
  roles <- rownames(M); np <- npop[roles]
  worst <- 0; pair <- "-"
  for (a in roles) for (b in roles) {
    if (a == b) next
    lhs <- np[a] * M[a, b]; rhs <- np[b] * M[b, a]
    rel <- abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-12)
    if (rel > worst) { worst <- rel; pair <- paste(a, b, sep = "<->") }
  }
  cat(sprintf("  %-15s max reciprocity asymmetry = %5.1f%%  (%s)\n",
              name, worst * 100, pair))
}
if (VERBOSE) {
  cat("\n── Reciprocity (N_i*M[i,j] vs N_j*M[j,i]; 0%% if balanced) ──\n")
  reciprocity_report(M_infant_room,  role_pop, "M_infant_room")
  reciprocity_report(M_toddler_room, role_pop, "M_toddler_room")
  reciprocity_report(M_household,    role_pop, "M_household")
}

# ── 11. SANITY CHECKS (M-appropriate) ─────────────────────────────────────────
# NOTE: unlike F, M entries are CONTACT COUNTS and routinely exceed 1.
#       The old ">1 is unexpected" check has been removed.
check_matrix <- function(M, name, soft_max = 50) {
  issues <- character(0)
  if (any(is.na(M)))  issues <- c(issues, "contains NA")
  if (any(M < 0))     issues <- c(issues, "contains negative values")
  if (any(M > soft_max))
    issues <- c(issues, sprintf("entry > %g contacts/day (check units)", soft_max))
  if (length(issues) == 0) {
    cat(sprintf("  OK  %-15s [range %.5f – %.5f]\n", name, min(M), max(M)))
  } else {
    warning(sprintf("  %s: %s", name, paste(issues, collapse = "; ")), call. = FALSE)
  }
}
cat("\n── Sanity checks ──\n")
check_matrix(M_infant_room,  "M_infant_room")
check_matrix(M_toddler_room, "M_toddler_room")
check_matrix(M_household,    "M_household")

# ── 12. EXPORT MATRICES AS CSV (OPTIONAL) ─────────────────────────────────────
write_matrix_csv <- function(M, filename) {
  df <- cbind(role = rownames(M), as.data.frame(M))
  write.csv(df, filename, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", filename))
}
# write_matrix_csv(M_infant_room,  "M_infant_room.csv")
# write_matrix_csv(M_toddler_room, "M_toddler_room.csv")
# write_matrix_csv(M_household,    "M_household.csv")

# ── 13. SUMMARY ───────────────────────────────────────────────────────────────
cat("\nDone. Three matrices available in the global environment:\n")
cat("  M_infant_room  [", paste(rownames(M_infant_room),  collapse = " x "), "]\n")
cat("  M_toddler_room [", paste(rownames(M_toddler_room), collapse = " x "), "]\n")
cat("  M_household    [", paste(rownames(M_household),    collapse = " x "), "]\n")
cat(sprintf("  Reciprocity balancing: %s\n\n",
            if (BALANCE_RECIPROCITY) "ON" else "OFF"))

# ── APPENDIX A: Role-to-band mapping ──────────────────────────────────────────
# Role     Age range   Epistorm-Mix 5yr bands            Band indices
# infant   0 – <1 yr   "0-4" (contact weight 1/5)        1
# toddler  1 –  4 yr   "0-4" (contact weight 4/5)        1
# staff    18 – 64 yr  "20-24" … "60-64"                 5:13
# elderly  65+         "65-69","70-74","75-79","80+"     14:17
# Excluded: "5-9","10-14","15-19" (school-age, not in model)
#
# ── APPENDIX B: Aggregation derivation (M-matrix) ─────────────────────────────
# M[i,j] = mean number of contacts a band-i participant has with band-j people.
# For super-groups I (participants) and J (contacts):
#   contacts of one band-i participant with all of J  = sum_{j in J} M[i,j]   (SUM)
#   a random I-participant is from band i w.p. N_i/N_I                        (rows)
#   => M[I,J] = sum_{i in I} (N_i/N_I) * sum_{j in J} M[i,j]
#   R: sum(pop_rows * rowSums(sub)) / sum(pop_rows)
# Population weighting is applied to ROWS ONLY; columns are an unweighted sum.
# This differs from an F-matrix, where the merged-contact entry is a per-capita
# probability and columns are NOT summed.
#
# Band-1 (infant/toddler) split:
#   contact (column): multiply by within-band population share (1/5 or 4/5)
#   participant (row): unchanged (band-average rate) -> infant & toddler rows
#                      are identical (see §4).
